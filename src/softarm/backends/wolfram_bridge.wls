(* Generic SymPy AST executor. Model-specific equations must never live here. *)

ClearAll[
  resetRequest, fromSymbol, fromGraph, toGraph, installDerivatives,
  requestParameter, requestTimeout, executeRequest, writeResponse, cseData
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
  optimized = TimeConstrained[
    Experimental`OptimizeExpression[expressions],
    60,
    $Failed
  ];
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
  {operation, expressions, variables, values, lower, upper, result,
   timeout, rows, columns, data},
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
        "optimize", FactorTerms /@ expressions,
        "differentiate",
          variables = fromGraph[request["parameters"]["variables"]];
          Flatten[D[expressions, {variables}]],
        "integrate",
          variables = {First[fromGraph[request["parameters"]["variable"]]]};
          lower = First[fromGraph[request["parameters"]["lower"]]];
          upper = First[fromGraph[request["parameters"]["upper"]]];
          Integrate[#, {First[variables], lower, upper},
            Assumptions -> And @@ assumptions] & /@ expressions,
        "substitute",
          variables = fromGraph[request["parameters"]["variables"]];
          values = fromGraph[request["parameters"]["values"]];
          (# /. Thread[variables -> values]) & /@ expressions,
        "rank",
          rows = requestParameter[request, "rows"];
          columns = requestParameter[request, "cols"];
          MatrixRank[Partition[expressions, columns]],
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
    "rank", <|"ok" -> True, "rank" -> result|>,
    "cse",
      If[AssociationQ[result],
        Join[<|"ok" -> True|>, result],
        <|"ok" -> False, "error" -> "unsupported optimized-expression shape"|>
      ],
    _, <|"ok" -> True, "expressions" -> toGraph[result]|>
  ]
];

writeResponse[response_Association] := (
  WriteString[$Output,
    ExportString[response, "RawJSON", "Compact" -> True], "\n"];
  Flush[$Output];
);

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
 If[Last[$CommandLine] === "--server",
  writeResponse[<|"ok" -> True, "version" -> $Version|>];
  standardInput = InputStream["stdin", 1];
  While[True,
    line = ReadLine[standardInput];
    If[line === EndOfFile, Break[]];
    request = Quiet[Check[ImportString[line, "RawJSON"], $Failed]];
    If[request === $Failed,
      writeResponse[<|"ok" -> False,
        "error" -> ("invalid JSON request: " <>
          StringTake[ToString[line, InputForm], UpTo[240]])|>];
      Continue[]
    ];
    If[Lookup[request, "operation", ""] === "close", Break[]];
    response = Quiet[Check[
      executeRequest[request],
      <|"ok" -> False, "error" -> "unhandled Wolfram evaluation error"|>
    ]];
    writeResponse[response];
  ];
  Exit[0]
];

(* Backward-compatible one-shot optimizer used by older callers. *)
args = Take[$CommandLine, -2];
If[Length[args] != 2,
  Print["expected input and output paths"];
  Exit[2]
];
payload = Import[args[[1]], "RawJSON"];
symbolByOriginal = <||>;
originalBySymbol = <||>;
fromLegacy[node_Association] := Module[{head = node["head"], values},
  Switch[head,
    "symbol",
      If[!KeyExistsQ[symbolByOriginal, node["name"]],
        AssociateTo[symbolByOriginal,
          node["name"] -> ("saSymbol" <> ToString[Length[symbolByOriginal] + 1])];
        AssociateTo[originalBySymbol,
          symbolByOriginal[node["name"]] -> node["name"]]
      ];
      Symbol[symbolByOriginal[node["name"]]],
    "integer", ToExpression[node["value"], InputForm],
    "rational",
      ToExpression[node["numerator"], InputForm]/
        ToExpression[node["denominator"], InputForm],
    "real", ToExpression[node["value"], InputForm],
    _,
      values = fromLegacy /@ Lookup[node, "args", {}];
      Switch[head,
        "add", Plus @@ values, "mul", Times @@ values,
        "pow", Power @@ values, "sin", Sin @@ values,
        "cos", Cos @@ values, "abs", Abs @@ values,
        _, Apply[Symbol[head], values]
      ]
  ]
];
toLegacy[value_Integer] :=
  <|"head" -> "integer", "value" -> ToString[value, InputForm]|>;
toLegacy[value_Rational] :=
  <|"head" -> "rational", "numerator" -> ToString[Numerator[value]],
    "denominator" -> ToString[Denominator[value]]|>;
toLegacy[value_Real] :=
  <|"head" -> "real", "value" -> ToString[value, InputForm]|>;
toLegacy[value_Symbol] := Module[{name = SymbolName[Unevaluated[value]]},
  <|"head" -> "symbol", "name" -> Lookup[originalBySymbol, name, name]|>
];
toLegacy[value_Plus] :=
  <|"head" -> "add", "args" -> (toLegacy /@ (List @@ value))|>;
toLegacy[value_Times] :=
  <|"head" -> "mul", "args" -> (toLegacy /@ (List @@ value))|>;
toLegacy[value_Power] :=
  <|"head" -> "pow", "args" -> (toLegacy /@ (List @@ value))|>;
toLegacy[value_Sin] := <|"head" -> "sin", "args" -> {toLegacy[value[[1]]]}|>;
toLegacy[value_Cos] := <|"head" -> "cos", "args" -> {toLegacy[value[[1]]]}|>;
toLegacy[value_Abs] := <|"head" -> "abs", "args" -> {toLegacy[value[[1]]]}|>;
toLegacy[value_] :=
  <|"head" -> SymbolName[Head[value]],
    "args" -> (toLegacy /@ (List @@ value))|>;
expressions = fromLegacy /@ payload["expressions"];
optimized = TimeConstrained[FactorTerms /@ expressions, 600, $Failed];
If[optimized === $Failed, Print["optimization timed out"]; Exit[3]];
Export[args[[2]],
  <|"expressions" -> (toLegacy /@ optimized)|>, "RawJSON"];
Exit[0];
];
