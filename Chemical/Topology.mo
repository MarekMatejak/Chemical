within Chemical;
package Topology
  extends Modelica.Icons.Package;

  model SplitterT1 "Splits a flow into two subflows"

    replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
      annotation (choices(
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Interfaces.Inlet inlet(redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-100,0})));
    Interfaces.Outlet outletA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,100})));
    Chemical.Interfaces.Outlet outletB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,-100})));
    SplitterN splitterN(redeclare package stateOfMatter = stateOfMatter,final N=2, final L=L)
      annotation (Placement(transformation(extent={{-28,-10},{-8,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(splitterN.inlet, inlet) annotation (Line(
        points={{-28,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet[2], outletA) annotation (Line(
        points={{-8,0.5},{-8,50},{0,50},{0,100}},
        color={158,66,200},
        thickness=0.5));
    connect(outletB, splitterN.outlet[1]) annotation (Line(
        points={{0,-100},{0,-0.5},{-8,-0.5}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,-100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{20,100},{60,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{20,-60},{60,-100}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SplitterT1;

  model SplitterT2 "Splits a flow into two subflows"

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Interfaces.Inlet inletProcess annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
    Interfaces.Outlet outletSubstance annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,100})));
    Interfaces.Outlet outletSubstance1 annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={100,0})));
    SplitterN splitterN(final N=2, final L=L)
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(splitterN.inlet, inletProcess) annotation (Line(
        points={{-40,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet[1], outletSubstance1) annotation (Line(
        points={{-20,-0.5},{-20,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(outletSubstance, splitterN.outlet[2]) annotation (Line(
        points={{0,100},{0,0.5},{-20,0.5}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{-20,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{60,-20},{100,-60}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SplitterT2;

  model JunctionT1 "2 to 1 T-Junction"

    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,0})));
    Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-100})));
    JunctionN junctionN(final N=2, final L=L,
       final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-34,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionN.outlet, outlet) annotation (Line(
        points={{-44,1.33227e-15},{-72,1.33227e-15},{-72,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(inletA, junctionN.inlets[1]) annotation (Line(
        points={{0,100},{0,0.5},{-24,0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionN.inlets[2], inletB) annotation (Line(
        points={{-24,-0.5},{0,-0.5},{0,-100}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,-100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{20,100},{60,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{20,-60},{60,-100}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionT1;

  model JunctionT2 "2 to 1 T-Junction"

    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,0})));
    Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={100,0})));
    JunctionN junctionN(final N=2, final L=L, final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-20,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionN.inlets[2], inletB) annotation (Line(
        points={{-10,-0.5},{36,-0.5},{36,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(inletA, junctionN.inlets[1]) annotation (Line(
        points={{0,100},{0,0.5},{-10,0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionN.outlet, outlet) annotation (Line(
        points={{-30,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{-20,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{80,-20},{120,-60}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionT2;

  model SplitterX "Splits a flow into three subflows"

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
    Chemical.Interfaces.Outlet outletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,100})));
    Chemical.Interfaces.Outlet outletB
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={3.55271e-15,-100})));
    Chemical.Interfaces.Outlet outletC annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={100,0})));
    SplitterN splitterN(final N=3, final L=L)
      annotation (Placement(transformation(extent={{-32,-10},{-12,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(splitterN.inlet, inlet) annotation (Line(
        points={{-32,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet[3], outletA) annotation (Line(
        points={{-12,0.666667},{0,0.666667},{0,100}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet[2], outletC) annotation (Line(
        points={{-12,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet[1], outletB)
      annotation (Line(
        points={{-12,-0.666667},{0,-0.666667},{0,-100},{3.55271e-15,-100}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,-100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-20,100},{-60,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{20,-60},{60,-100}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{50, 60},{90, 20}},
            textColor={175,175,175},
            textString="C")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SplitterX;

  model JunctionX1 "2 to 2 X-Junction"

    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Outlet outleta annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,0})));
    Chemical.Interfaces.Outlet outletb annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,-100})));
    Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={100,0})));
    JunctionNM junctionNM(N=2, M=2, final L=L, final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={0,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(inletB, junctionNM.inlets[2]) annotation (Line(
        points={{100,0},{54,0},{54,-0.5},{10,-0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionNM.outlets[1], outleta) annotation (Line(
        points={{-10,0.5},{-54,0.5},{-54,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(inletA, junctionNM.inlets[1]) annotation (Line(
        points={{0,100},{0,40},{40,40},{40,0.5},{10,0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(outletb, junctionNM.outlets[2]) annotation (Line(
        points={{0,-100},{0,-40},{-40,-40},{-40,-0.5},{-10,-0.5}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,-100},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{-20,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{50,20},{90,60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{-60,-20},{-100,-60}},
            textColor={175,175,175},
            textString="a"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="b")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Medium package used in the component. Make sure it is the same one as all the components connected to all fluid ports are using. </p>
</html>"));
  end JunctionX1;

  model JunctionX2 "2 to 2 X-Junction"

    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Outlet outleta annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,0})));
    Chemical.Interfaces.Outlet outletb annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={100,0})));
    Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-100})));
    JunctionNM junctionNM(N=2, M=2, final L=L, final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={0,20})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionNM.outlets[1], outleta) annotation (Line(
        points={{-0.5,10},{-1,10},{-1,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionNM.outlets[2], outletb) annotation (Line(
        points={{0.5,10},{1,10},{1,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionNM.inlets[1], inletA) annotation (Line(
        points={{-0.5,30},{-0.5,66},{0,66},{0,100}},
        color={158,66,200},
        thickness=0.5));
    connect(inletB, junctionNM.inlets[2]) annotation (Line(
        points={{0,-100},{0,-40},{40,-40},{40,52},{0.5,52},{0.5,30}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,-100},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{-20,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{-60,-20},{-100,-60}},
            textColor={175,175,175},
            textString="a"),
          Text(
            extent={{50,20},{90,60}},
            textColor={175,175,175},
            textString="b")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionX2;

  model JunctionX3 "3 to 1 X-Junction"

    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={-100,0})));
    Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=-90,
          origin={0,100})));
    Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={100,0})));
    Chemical.Interfaces.Inlet inletC annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-100})));
    JunctionN junctionN(final N=3, final L=L,
      final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-40,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionN.inlets[2], inletB) annotation (Line(
        points={{-30,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(inletC, junctionN.inlets[3]) annotation (Line(
        points={{0,-100},{0,-0.666667},{-30,-0.666667}},
        color={158,66,200},
        thickness=0.5));
    connect(inletA, junctionN.inlets[1]) annotation (Line(
        points={{0,100},{0,0.666667},{-30,0.666667}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionN.outlet, outlet) annotation (Line(
        points={{-50,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,100}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{0,-100}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{-20,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{50,20},{90,60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="C")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionX3;

  model SplitterN "Splitter with one inlet and N outlets"

    replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
      annotation (choices(
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    parameter Integer N(min=1) = 1 "Number of outputs";
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Interfaces.Inlet inlet(redeclare package stateOfMatter = stateOfMatter) "inlet" annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));
    Interfaces.Outlet outlet[N](redeclare package stateOfMatter = stateOfMatter) "vector of N outlets"
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential r_mix;

  equation
    der(inlet.n_flow)*L = inlet.r - r_mix;

    for i in 1:N loop
      der(outlet[i].n_flow)*L = outlet[i].r - r_mix;
      outlet[i].definition = inlet.definition;
      outlet[i].solution = inlet.solution;

      outlet[i].state.u = inlet.state.u;
      outlet[i].state.h = inlet.state.h;
    end for;

    sum(outlet.n_flow) + inlet.n_flow = 0;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{0,0},{96,10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{96,-10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{-100,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{90,80},{50,40}},
            textColor={175,175,175},
            textString="%N")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SplitterN;

  model JunctionN "Junction with N inlets and one outlet"

    parameter Integer N(min=1) = 1 "Number of inlets";
    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlets[N] "vector of N inlets"
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));
    Chemical.Interfaces.Outlet outlet "outlet"
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

    // these are needed by DynamicJunctionN
    output Real w[N](each unit="1") "regularized weighting factor for specific enthalpy";

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential u[N]=inlets.state.u "(steady molar-flow) electro-chemical potential";
    Modelica.Units.SI.MolarEnthalpy h[N]=inlets.state.h "molar enthapy at inlets";

    Modelica.Units.SI.ChemicalPotential u_mix "(steady mass-flow) electro-chemical potential at the outlet";

    Modelica.Units.SI.ChemicalPotential r_mix "inertial electro-chemical potential divided by R*T at outlet";
    Modelica.Units.SI.MolarEnthalpy h_mix "molar enthalpy at outlet";

    Modelica.Units.SI.ChemicalPotential r_in[N];

  equation
    sum(inlets.n_flow) + outlet.n_flow = 0;

    for i in 1:N loop
      der(inlets[i].n_flow) * L = inlets[i].r - r_in[i];

      u[i] + r_in[i] = u_mix + r_mix;
      w[i] = (abs(inlets[i].n_flow)+n_flow_eps) / (sum(abs(inlets.n_flow))+N*n_flow_eps);

      //assert(outlet.definition.DGH >= inlets[i].definition,"Topology connections allowed only for the same substance.");
      //assert(outlet.solution.T == inlets[i].solution.T,"Topology connections allowed only for the same chemical solution.");
    end for;
    der(outlet.n_flow) * L =  outlet.r - r_mix;

    u_mix = sum(w.*u);
    h_mix = sum(w.*h);

    outlet.state.u = u_mix;
    outlet.state.h = h_mix;

    outlet.definition = inlets[1].definition;
    outlet.solution = inlets[1].solution;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{-100,10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{-100,-10}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-90,80},{-50,40}},
            textColor={175,175,175},
            textString="%N")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionN;

  model JunctionNM "Junction with N inlets and M outlets"

    parameter Integer N(min=1) = 1 "Number of inputs";
    parameter Integer M(min=1) = 1 "Number of outputs";
    parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlets[N] "vector of N inlets" annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-100,0}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-100,0})));
    Chemical.Interfaces.Outlet outlets[M] "vector of N outlets" annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={100,0}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={100,0})));
    SplitterN splitterN(
      final N=M,
      final L=L)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={18,0})));
    JunctionN junctionN(
      final N=N,
      final L=L,
      final n_flow_eps=n_flow_eps)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-14,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionN.inlets, inlets) annotation (Line(
        points={{-24,0},{-62,0},{-62,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionN.outlet, splitterN.inlet) annotation (Line(
        points={{-4,0},{8,0}},
        color={158,66,200},
        thickness=0.5));
    connect(splitterN.outlet, outlets) annotation (Line(
        points={{28,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{-100,10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{-100,-10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{96,10}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{0,0},{96,-10}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-90,80},{-50,40}},
            textColor={175,175,175},
            textString="%N"),
          Text(
            extent={{90,80},{50,40}},
            textColor={175,175,175},
            textString="%M")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end JunctionNM;
annotation (Documentation(revisions="<html>
<p><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>

</html>"), Icon(graphics={
        Line(
          points={{-80,0},{12,0}},
          color={158,66,200},
          thickness=0.5),
        Line(
          points={{12,0},{80,-80}},
          color={158,66,200},
          thickness=0.5),
        Line(
          points={{12,0},{80,80}},
          color={158,66,200},
          thickness=0.5),
        Ellipse(
          extent={{6,6},{18,-6}},
          lineColor={158,66,200},
          fillColor={194,138,221},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5)}));
end Topology;
