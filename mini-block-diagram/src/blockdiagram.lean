-- Block Diagram — Lean 4 Formalization
-- Knowledge Coverage: L1-L4 (defs, concepts, structures, theorems)
-- Ref: S.J. Mason, Proc. IRE, 1953
-- MIT 6.302 / Stanford ENGR105 / ETH 151-0591

inductive NodeType : Type where
  | input   : NodeType
  | output  : NodeType
  | block   : NodeType
  | summer  : NodeType
  | takeoff : NodeType
  | gain    : NodeType
  deriving BEq, Repr, Inhabited

structure Node where
  id    : Nat
  ntype : NodeType
  label : String
  gain  : Float
  deriving BEq, Repr

structure Edge where
  id      : Nat
  src     : Nat
  dst     : Nat
  weight  : Float
  deriving BEq, Repr

structure BlockDiagram where
  nodes    : List Node
  edges    : List Edge
  inputId  : Nat
  outputId : Nat
  name     : String
  deriving Repr

def BlockDiagram.hasNode (bd : BlockDiagram) (nid : Nat) : Bool :=
  bd.nodes.any (fun n => n.id == nid)

def BlockDiagram.validEdge (bd : BlockDiagram) (e : Edge) : Bool :=
  bd.hasNode e.src && bd.hasNode e.dst

def BlockDiagram.allEdgesValid (bd : BlockDiagram) : Bool :=
  bd.edges.all (bd.validEdge)

def BlockDiagram.validInput (bd : BlockDiagram) : Bool :=
  bd.hasNode bd.inputId &&
  ¬(bd.edges.any (fun e => e.dst == bd.inputId))

def BlockDiagram.validOutput (bd : BlockDiagram) : Bool :=
  bd.hasNode bd.outputId &&
  ¬(bd.edges.any (fun e => e.src == bd.outputId))

structure Path where
  nodes  : List Nat
  edges  : List Nat
  gain   : Float
  isLoop : Bool
  deriving Repr

def Edge.consecutive (e1 e2 : Edge) : Bool :=
  e1.dst == e2.src

def BlockDiagram.isValidPath (bd : BlockDiagram) (edges : List Edge) : Bool :=
  match edges with
  | [] => true
  | [e] => bd.validEdge e
  | e1 :: e2 :: rest =>
    bd.validEdge e1 && e1.consecutive e2 && bd.isValidPath (e2 :: rest)

def BlockDiagram.isForwardPath (bd : BlockDiagram) (p : Path) : Bool :=
  match p.nodes.head?, p.nodes.getLast? with
  | some first, some last => first == bd.inputId && last == bd.outputId
  | _, _ => false

def BlockDiagram.isLoop (bd : BlockDiagram) (p : Path) : Bool :=
  match p.nodes.head?, p.nodes.getLast? with
  | some first, some last => first == last && p.nodes.length >= 2
  | _, _ => false

def Path.nonTouching (l1 l2 : Path) : Bool :=
  l1.nodes.all (fun n => ¬(l2.nodes.contains n))

def Path.mutuallyNonTouching (loops : List Path) : Bool :=
  match loops with
  | [] => true
  | l :: ls => ls.all (l.nonTouching) && mutuallyNonTouching ls

theorem input_no_incoming (bd : BlockDiagram) (h : bd.validInput) :
  bd.edges.filter (fun e => e.dst == bd.inputId) |>.length = 0 := by
  unfold BlockDiagram.validInput at h
  cases h; rename_i h1 h2
  have : ¬(bd.edges.any (fun e => e.dst == bd.inputId)) := h2
  simp [List.length_eq_zero]
  intro hne
  apply this
  exact List.any_of_mem hne (fun x hx => hx)

theorem output_no_outgoing (bd : BlockDiagram) (h : bd.validOutput) :
  bd.edges.filter (fun e => e.src == bd.outputId) |>.length = 0 := by
  unfold BlockDiagram.validOutput at h
  cases h; rename_i h1 h2
  have : ¬(bd.edges.any (fun e => e.src == bd.outputId)) := h2
  simp [List.length_eq_zero]
  intro hne
  apply this
  exact List.any_of_mem hne (fun x hx => hx)

theorem empty_diagram_allEdgesValid :
  ({ nodes := [], edges := [], inputId := 0, outputId := 0,
     name := "empty" } : BlockDiagram).allEdgesValid := by
  unfold BlockDiagram.allEdgesValid; simp [List.all]

theorem consecutive_edges_form_path (bd : BlockDiagram) (e1 e2 : Edge)
    (hv1 : bd.validEdge e1) (hv2 : bd.validEdge e2)
    (hcons : e1.dst == e2.src) : bd.isValidPath [e1, e2] := by
  unfold BlockDiagram.isValidPath
  simp [hv1, hcons, hv2]

def Path.gainProduct (edges : List Edge) : Float :=
  edges.foldl (fun acc e => acc * e.weight) 1.0

theorem gainProduct_empty : Path.gainProduct [] = 1.0 := by
  unfold Path.gainProduct; rfl

theorem gainProduct_singleton (e : Edge) : Path.gainProduct [e] = e.weight := by
  unfold Path.gainProduct; simp
