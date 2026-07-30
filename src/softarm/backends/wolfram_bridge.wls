(* math.exe 15 on Windows leaves $ScriptCommandLine empty for -script;
   the two explicit paths are always the final command-line arguments. *)
args = Take[$CommandLine, -2];
If[Length[args] != 2, Print["expected input and output paths"]; Exit[2]];

symbolByOriginal = <||>;
originalBySymbol = <||>;
fromSymbol[name_String] := Module[{safe},
  If[KeyExistsQ[symbolByOriginal, name], Return[Symbol[symbolByOriginal[name]]]];
  safe = "saSymbol" <> ToString[Length[symbolByOriginal] + 1];
  AssociateTo[symbolByOriginal, name -> safe];
  AssociateTo[originalBySymbol, safe -> name];
  Symbol[safe]
];

fromNode[node_Association] := Module[{head = node["head"], values},
  Switch[head,
    "symbol", fromSymbol[node["name"]],
    "integer", ToExpression[node["value"], InputForm],
    "rational", ToExpression[node["numerator"], InputForm]/ToExpression[node["denominator"], InputForm],
    "real", ToExpression[node["value"], InputForm],
    _,
      values = fromNode /@ Lookup[node, "args", {}];
      Switch[head,
        "add", Plus @@ values,
        "mul", Times @@ values,
        "pow", Power @@ values,
        "sin", Sin @@ values,
        "cos", Cos @@ values,
        "abs", Abs @@ values,
        _, Apply[Symbol[head], values]
      ]
  ]
];

toNode[value_Integer] := <|"head" -> "integer", "value" -> ToString[value, InputForm]|>;
toNode[value_Rational] := <|"head" -> "rational", "numerator" -> ToString[Numerator[value]], "denominator" -> ToString[Denominator[value]]|>;
toNode[value_Real] := <|"head" -> "real", "value" -> ToString[value, InputForm]|>;
toNode[value_Symbol] := Module[{name = SymbolName[Unevaluated[value]]},
  <|"head" -> "symbol", "name" -> Lookup[originalBySymbol, name, name]|>
];
toNode[value_Plus] := <|"head" -> "add", "args" -> (toNode /@ (List @@ value))|>;
toNode[value_Times] := <|"head" -> "mul", "args" -> (toNode /@ (List @@ value))|>;
toNode[value_Power] := <|"head" -> "pow", "args" -> (toNode /@ (List @@ value))|>;
toNode[value_Sin] := <|"head" -> "sin", "args" -> {toNode[value[[1]]]}|>;
toNode[value_Cos] := <|"head" -> "cos", "args" -> {toNode[value[[1]]]}|>;
toNode[value_Abs] := <|"head" -> "abs", "args" -> {toNode[value[[1]]]}|>;
toNode[value_] := <|"head" -> SymbolName[Head[value]], "args" -> (toNode /@ (List @@ value))|>;

payload = Import[args[[1]], "RawJSON"];
expressions = fromNode /@ payload["expressions"];
optimized = TimeConstrained[FactorTerms /@ expressions, 600, expressions];
Export[args[[2]], <|"expressions" -> (toNode /@ optimized)|>, "RawJSON"];
Exit[0];
