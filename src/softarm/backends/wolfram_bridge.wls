(* Generic SymPy AST executor. Model-specific equations must never live here. *)

ClearAll[
  resetRequest, fromSymbol, fromGraph, toGraph, installDerivatives,
  requestParameter, requestTimeout, executeRequest, cseData
];

resetRequest[] := (
  symbolByOriginal = <||>;
  originalBySymbol = <||>;
  assumptions = {};
);

fromSymbol[node_Association] := Module[
  {name = node["name"], safe, properties},
  If[KeyExistsQ[symbolByOriginal, name],
    Return[Symbol[symbolByOriginal[name]]]
  ];
  safe = "saSymbol" <> ToString[Length[symbolByOriginal] + 1];
  AssociateTo[symbolByOriginal, name -> safe];
  AssociateTo[originalBySymbol, safe -> name];
  properties = Lookup[node, "assumptions", <||>];
  If[TrueQ[Lookup[properties, "positive", False]],
    AppendTo[assumptions, Symbol[safe] > 0],
    If[TrueQ[Lookup[properties, "nonnegative", False]],
      AppendTo[assumptions, Symbol[safe] >= 0],
      If[TrueQ[Lookup[properties, "real", False]],
        AppendTo[assumptions, Element[Symbol[safe], Reals]]
      ]
    ]
  ];
  Symbol[safe]
];

fromGraph[graph_Association] := Module[{values = {}, node, args, value},
  Do[
    node = item;
    args = (values[[# + 1]] &) /@ Lookup[node, "args", {}];
    value = Switch[node["head"],
      "symbol", fromSymbol[node],
      "integer", ToExpression[node["value"], InputForm],
      "rational",
        ToExpression[node["numerator"], InputForm]/
          ToExpression[node["denominator"], InputForm],
      "real", ToExpression[node["value"], InputForm],
      "add", Plus @@ args,
      "mul", Times @@ args,
      "pow", Power @@ args,
      "sin", Sin @@ args,
      "cos", Cos @@ args,
      "abs", Abs @@ args,
      _, Apply[Symbol[node["head"]], args]
    ];
    AppendTo[values, value],
    {item, graph["nodes"]}
  ];
  (values[[# + 1]] &) /@ graph["roots"]
];

toGraph[expressions_List] := Module[
  {nodes = {}, indices = <||>, visit, roots},
  visit[value_] := Module[{key, cached, args, name, node, index},
    key = ToString[HoldComplete[value], InputForm];
    cached = Lookup[indices, key, Missing["NotFound"]];
    If[!MissingQ[cached], Return[cached]];
    args = visit /@ (List @@ value);
    node = Which[
      IntegerQ[value],
        <|"head" -> "integer", "value" -> ToString[value, InputForm]|>,
      Head[value] === Rational,
        <|"head" -> "rational",
          "numerator" -> ToString[Numerator[value]],
          "denominator" -> ToString[Denominator[value]]|>,
      Head[value] === Real,
        <|"head" -> "real", "value" -> ToString[value, InputForm]|>,
      Head[value] === Symbol,
        name = SymbolName[Unevaluated[value]];
        <|"head" -> "symbol",
          "name" -> Lookup[originalBySymbol, name, name]|>,
      Head[value] === Plus, <|"head" -> "add", "args" -> args|>,
      Head[value] === Times, <|"head" -> "mul", "args" -> args|>,
      Head[value] === Power, <|"head" -> "pow", "args" -> args|>,
      Head[value] === Sin, <|"head" -> "sin", "args" -> args|>,
      Head[value] === Cos, <|"head" -> "cos", "args" -> args|>,
      Head[value] === Abs, <|"head" -> "abs", "args" -> args|>,
      True, <|"head" -> SymbolName[Head[value]], "args" -> args|>
    ];
    index = Length[nodes];
    AppendTo[nodes, node];
    AssociateTo[indices, key -> index];
    index
  ];
  roots = visit /@ expressions;
  <|"nodes" -> nodes, "roots" -> roots|>
];

installDerivatives[rules_List] := Scan[
  Function[rule,
    With[
      {head = Symbol[rule["head"]],
       derivative = Symbol[rule["derivative"]]},
      head /: Derivative[1][head] = derivative
    ]
  ],
  rules
];

requestParameter[request_, name_] := Module[{value},
  value = Lookup[Lookup[request, "parameters", <||>], name, Missing[name]];
  If[MissingQ[value], Return[value]];
  If[AssociationQ[value] && KeyExistsQ[value, "nodes"],
    First[fromGraph[value]],
    value
  ]
];

requestTimeout[request_] := N[Lookup[request, "timeout", 600.0]];

cseData[expressions_List, prefix_String] := Module[
  {optimized, held, tree, parsed, locals, body, assignments, reduced,
   replacementRules, names},
  optimized = Quiet[Check[
    Experimental`OptimizeExpression[expressions],
    $Failed
  ]];
  If[optimized === $Failed, Return[$Failed]];
  held = optimized /.
    Experimental`OptimizedExpression[value_] :> HoldComplete[value];
  If[Head[held] =!= HoldComplete, Return[$Failed]];
  tree = held /. {Block -> List, CompoundExpression -> List, Set -> Rule};
  parsed = Quiet[Check[ReleaseHold[tree], $Failed]];
  If[parsed === $Failed || !ListQ[parsed] || Length[parsed] != 2,
    Return[$Failed]
  ];
  locals = parsed[[1]];
  body = parsed[[2]];
  If[!ListQ[locals], Return[$Failed]];
  If[!ListQ[body], body = {body}];
  assignments = Select[Most[body], MatchQ[#, _Rule] &];
  reduced = Last[body];
  If[Length[assignments] != Length[locals] || !ListQ[reduced],
    Return[$Failed]
  ];
  names = Table[prefix <> ToString[index - 1], {index, Length[locals]}];
  replacementRules = Thread[locals -> (Symbol /@ names)];
  <|
    "replacement_names" -> names,
    "replacements" -> toGraph[(Last /@ assignments) /. replacementRules],
    "expressions" -> toGraph[reduced /. replacementRules]
  |>
];

executeRequest[request_Association] := Module[
  {operation, expressions, variables, result, timeout, data},
  If[Lookup[request, "schema", 0] != 2,
    Return[<|"ok" -> False, "error" -> "unsupported protocol schema"|>]
  ];
  resetRequest[];
  expressions = fromGraph[request["expressions"]];
  installDerivatives[Lookup[request, "derivatives", {}]];
  operation = request["operation"];
  timeout = requestTimeout[request];
  result = TimeConstrained[
    Assuming[
      And @@ assumptions,
      Switch[operation,
        "factor_terms", FactorTerms /@ expressions,
        "differentiate",
          variables = fromGraph[request["parameters"]["variables"]];
          Flatten[D[expressions, {variables}]],
        "cse",
          data = cseData[
            expressions,
            ToString[requestParameter[request, "prefix"]]
          ];
          data,
        _, $Failed
      ]
    ],
    timeout,
    $Failed
  ];
  If[result === $Failed,
    Return[<|"ok" -> False,
      "error" -> ("operation failed or timed out: " <> operation)|>]
  ];
  Switch[operation,
    "cse",
      If[AssociationQ[result],
        Join[<|"ok" -> True|>, result],
        <|"ok" -> False, "error" -> "unsupported optimized-expression shape"|>
      ],
    _, <|"ok" -> True, "expressions" -> toGraph[result]|>
  ]
];

softarmExecuteJSON[text_String] := Module[{request},
  request = Quiet[Check[ImportString[text, "RawJSON"], $Failed]];
  If[request === $Failed,
    Return[ExportString[
      <|"ok" -> False, "error" -> "invalid JSON request"|>,
      "RawJSON", "Compact" -> True
    ]]
  ];
  ExportString[executeRequest[request], "RawJSON", "Compact" -> True]
];

If[!TrueQ[$SoftArmLibraryMode],
  Print["wolfram_bridge.wls is a library loaded by the Python backend"];
  Exit[2]
];
