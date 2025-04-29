within Chemical;
package Topology "Junctions and Connectors for undirected chemical simulation"
  extends Modelica.Icons.Package;

  model JunctionMN "Generalized junction/splitter for undirected flow"


    parameter Integer N(min=0) = 1 "Number of rears";
    parameter Integer M(min=0) = 1 "Number of fors";

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance for each branch" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of molar flow rate"
      annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rears[N] "Rear ports"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore fores[M] "Fore ports"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={100,0})));

    Modelica.Units.SI.ChemicalPotential u_mix "mixing u assuming positive molarFlow";

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential rs[M + N] "rs of the connectors internal of L";
    Modelica.Units.SI.ChemicalPotential r_mix "mixed r assuming positive molarflow";

    Modelica.Units.SI.MolarFlowRate inflows[M + N] "molarFlows assuming positive molarFlow";

    Modelica.Units.SI.ChemicalPotential us[M + N] "potential of medium entering";
    Modelica.Units.SI.ChemicalPotential us_out[M + N] "potential of medium exeting";

    Modelica.Units.SI.MolarEnthalpy hs[M + N] "molar enthalpy of medium entering";
    Modelica.Units.SI.MolarEnthalpy hs_out[M + N] "molar enthalpy of medium exeting";

  equation
    // rears are 1:N
    der(rears.n_flow)*L = rears.r-rs[1:N];
    for i in 1:N loop
       // inputs ports
      us[i] = rears[i].state_forwards.u;
      hs[i] = rears[i].state_forwards.h;
      // inflows, u_equation and set output
      inflows[i] = max(rears[i].n_flow, n_flow_reg);
      Chemical.Utilities.Internal.regStep(
          rears[i].n_flow,
          us[i],
          us_out[i],
          n_flow_reg) + rs[i] = u_mix + r_mix;
      rears[i].state_rearwards = Chemical.Interfaces.SubstanceState(u=us_out[i],h=hs_out[i]);

      rears[i].solution_rearwards = fores[1].solution_rearwards;
    end for;

    // fores are N+1:end
    der(fores.n_flow)*L = fores.r-rs[N+1:end];
    for i in 1:M loop
       // inputs ports
      us[N+i] = fores[i].state_rearwards.u;
      hs[N+i] = fores[i].state_rearwards.h;
      // inflows, u_equation and set output
      inflows[N+i] = max(fores[i].n_flow, n_flow_reg);
      Chemical.Utilities.Internal.regStep(
          fores[i].n_flow,
          us[N + i],
          us_out[N + i],
          n_flow_reg) + rs[N + i] = u_mix + r_mix;
      fores[i].state_forwards = Chemical.Interfaces.SubstanceState(u=us_out[N+i],h=hs_out[N+i]);

      fores[i].definition = rears[1].definition;
      fores[i].solution_forwards = rears[1].solution_forwards;

    end for;

    // Mass balance
    sum(rears.n_flow) + sum(fores.n_flow) = 0;
    //compute u_mix for rs computation
    u_mix =(us*inflows)/sum(inflows);

    // compute output quantities
    for i in 1:M+N loop
      hs_out[i] = (hs*inflows - hs[i]*inflows[i]) /(sum(inflows) - inflows[i]);
      us_out[i] = (us*inflows - us[i]*inflows[i]) /(sum(inflows) - inflows[i]);
    end for;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-84,0},{84,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{60,20},{100,60}},
            textColor={175,175,175},
            textString="%M"),
          Text(
            extent={{-100,20},{-60,60}},
            textColor={175,175,175},
            textString="%N")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This model represents a generalized junction/splitter for undirected flow with N rear and M fore ports. </u>
<u>Note that in the undirected case a distinction between junction and splitter is not possible, since the flow direction is unknown in advance. </u>
</html>"));
  end JunctionMN;

  model ConnectForeFore "Undirected connector with fore and fore"
    extends Chemical.Topology.Internal.PartialSubstanceAndSolutionSource;
   /* */
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Fore fore_a
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.Fore fore_b
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    fore_b.state_forwards = fore_a.state_rearwards;
    fore_a.state_forwards = fore_b.state_rearwards;
    fore_a.n_flow + fore_b.n_flow = 0;
    L*der(fore_a.n_flow) = fore_a.r - fore_b.r;

    fore_b.definition=substanceData;
    fore_a.definition=substanceData;
    fore_b.solution_forwards=solutionState;
    fore_a.solution_forwards=solutionState;

    annotation (Icon(
        graphics={
          Line(
            points={{-20,0},{20,0}},
            color={158,66,200},
            thickness=0.5)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect two undirected rear ports. </u>
<u>Basically the connector switches the names of output and input of the two ports.</u>
</html>"));
  end ConnectForeFore;

  model ConnectRearRear "Undirected connector with rear and rear"


    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rear_a
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.Rear rear_b
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    rear_b.state_rearwards = rear_a.state_forwards;
    rear_a.state_rearwards = rear_b.state_forwards;
    rear_a.n_flow + rear_b.n_flow = 0;
    L*der(rear_a.n_flow) = rear_a.r - rear_b.r;

    rear_b.solution_rearwards = rear_a.solution_forwards;
    rear_a.solution_rearwards = rear_b.solution_forwards;
    annotation (Icon(
        graphics={
          Line(
            points={{-20,0},{20,0}},
            color={158,66,200},
            thickness=0.5)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect two undirected fore ports. </u>
<u>Basically the connector switches the names of output and input of the two ports.</u>
</html>"));
  end ConnectRearRear;

  model JunctionRFF "Junction with rear and two fores"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rear
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore foreA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Fore foreB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=2,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionMN.rears[1], rear) annotation (Line(
        points={{-50,0},{-76,0},{-76,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[1], foreB) annotation (Line(
        points={{-30,-0.5},{0,-0.5},{0,-100}},
        color={158,66,200},
        thickness=0.5));
    connect(foreA, junctionMN.fores[2]) annotation (Line(
        points={{0,100},{0,0.5},{-30,0.5}},
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
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>
Junction with a rear and two fores in a lying T shape.
</u>
</html>"));
  end JunctionRFF;

  model JunctionRFF2 "Junction with rear and two fores"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rear
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore foreA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Fore foreB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=2,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionMN.rears[1], rear) annotation (Line(
        points={{-50,0},{-76,0},{-76,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[2], foreA) annotation (Line(
        points={{-30,0.5},{-30,1},{0,1},{0,100}},
        color={158,66,200},
        thickness=0.5));
    connect(foreB, junctionMN.fores[1]) annotation (Line(
        points={{100,0},{20,0},{20,-0.5},{-30,-0.5}},
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
            extent={{20,100},{60,60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{60,-20},{100,-60}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with a rear and two fores in a standing T shape.</u>
</html>"));
  end JunctionRFF2;

  model JunctionRRF "Junction with two rears and a fore"

         parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Fore fore
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.Rear rearA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Rear rearB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=1,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{30,-10},{50,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(rearB, junctionMN.rears[1]) annotation (Line(
        points={{0,-100},{0,-2},{30,-2},{30,-0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(rearA, junctionMN.rears[2]) annotation (Line(
        points={{0,100},{0,0.5},{30,0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[1], fore) annotation (Line(
        points={{50,0},{76,0},{76,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{100,0},{0,0}},
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
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with two rears and a fore in a lying T shape.</u>
</html>"));
  end JunctionRRF;

  model JunctionRRF2 "Junction with two rears and a fore"

         parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    JunctionMN junctionMN(
      final M=1,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{30,-10},{50,10}})));
    Chemical.Interfaces.Rear rearA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    Chemical.Interfaces.Rear rearB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore fore
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(junctionMN.rears[1], rearA) annotation (Line(
        points={{30,-0.5},{0,-0.5},{0,-100}},
        color={158,66,200},
        thickness=0.5));
    connect(fore, junctionMN.fores[1]) annotation (Line(
        points={{100,0},{76,0},{76,0},{50,0}},
        color={158,66,200},
        thickness=0.5));
    connect(rearB, junctionMN.rears[2]) annotation (Line(
        points={{-100,0},{-36,0},{-36,0.5},{30,0.5}},
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
            extent={{20,-100},{60,-60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{60,20},{100,60}},
            textColor={175,175,175},
            textString="B")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with two rears and a fore in a standing T shape.</u>
</html>"));
  end JunctionRRF2;

  model JunctionRFFF "Junction with a rear and three fores"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rear
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore foreA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Fore foreB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.Fore foreC
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=3,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(junctionMN.rears[1], rear) annotation (Line(
        points={{-50,0},{-76,0},{-76,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[3], foreA) annotation (Line(
        points={{-30,0.666667},{0,0.666667},{0,100}},
        color={158,66,200},
        thickness=0.5));
    connect(foreB, junctionMN.fores[2]) annotation (Line(
        points={{100,0},{-30,0}},
        color={158,66,200},
        thickness=0.5));
    connect(foreC, junctionMN.fores[1]) annotation (Line(
        points={{0,-100},{0,-0.666667},{-30,-0.666667}},
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
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="C")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with a rear and three fores in a x shape.</u>
</html>"));
  end JunctionRFFF;

  model JunctionRRFF "Junction with two rears and two fores"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear reara
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Rear rearb
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Fore foreB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.Fore foreA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=2,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(foreB, junctionMN.fores[2]) annotation (Line(
        points={{100,0},{40,0},{40,0.5},{10,0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(foreA, junctionMN.fores[1]) annotation (Line(
        points={{0,-100},{0,-20},{10,-20},{10,-0.5}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.rears[1], reara) annotation (Line(
        points={{-10,-0.5},{-40,-0.5},{-40,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(rearb, junctionMN.rears[2]) annotation (Line(
        points={{0,100},{0,20},{-10,20},{-10,0.5}},
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
            textString="b"),
          Text(
            extent={{50,20},{90,60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="A"),
          Text(
            extent={{-50,-20},{-90,-60}},
            textColor={175,175,175},
            textString="a")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with two rears and two fores in a x shape.</u>
</html>"));
  end JunctionRRFF;

  model JunctionRRFF2 "Junction with two rears and two fores"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear reara
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Rear rearb
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    Chemical.Interfaces.Fore foreA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Fore foreB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=2,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{10,-10},{-10,10}},
          rotation=90,
          origin={0,10})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(rearb, junctionMN.rears[1]) annotation (Line(
        points={{0,-100},{0,-20},{20,-20},{20,20},{0.5,20}},
        color={158,66,200},
        thickness=0.5));
    connect(reara, junctionMN.rears[2]) annotation (Line(
        points={{0,100},{0,20},{-0.5,20}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[2], foreA) annotation (Line(
        points={{-0.5,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(foreB, junctionMN.fores[1]) annotation (Line(
        points={{100,0},{0.5,0}},
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
            textString="a"),
          Text(
            extent={{50,20},{90,60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="b"),
          Text(
            extent={{-50,-20},{-90,-60}},
            textColor={175,175,175},
            textString="A")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with two rears and two fores in a x shape.</u>
</html>"));
  end JunctionRRFF2;

  model JunctionRRRF "Junction with tree rears and a fore"

        parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Rear rearA
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.Fore fore
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=1,
      final N=3,
      final n_flow_reg=n_flow_reg,
      final L=L)
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Chemical.Interfaces.Rear rearB
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.Rear rearC
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(rearA, junctionMN.rears[2]) annotation (Line(
        points={{0,100},{0,0}},
        color={158,66,200},
        thickness=0.5));
    connect(rearB, junctionMN.rears[1]) annotation (Line(
        points={{-100,0},{0,0},{0,-0.666667}},
        color={158,66,200},
        thickness=0.5));
    connect(rearC, junctionMN.rears[3]) annotation (Line(
        points={{0,-100},{0,0.666667}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionMN.fores[1], fore) annotation (Line(
        points={{20,0},{60,0},{60,0},{100,0}},
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
            extent={{-50,-20},{-90,-60}},
            textColor={175,175,175},
            textString="B"),
          Text(
            extent={{60,-100},{20,-60}},
            textColor={175,175,175},
            textString="C")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Junction with three rears and a fore in a x shape.</u>
</html>"));
  end JunctionRRRF;

  package Internal
    partial model PartialSubstanceAndSolutionSource
      "Substance properties for components, where the substance is connected with the solution"
      import Chemical;


      parameter Chemical.Interfaces.Definition substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable=not useRear));

      parameter Boolean useSolution = false "Use solution connector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Chemical.Interfaces.SolutionState solutionParam = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible)
      "Constant chemical solution state if not from rear or input"
        annotation (Dialog(enable=not useSolution and not useRear));

      Chemical.Interfaces.SolutionPort solution(
          T=solutionState.T,
          p=solutionState.p,
          v=solutionState.v,
          n=solutionState.n,
          m=solutionState.m,
          V=solutionState.V,
          G=solutionState.G,
          Q=solutionState.Q,
          I=solutionState.I,
          i=0,
          dH=0,
          dV=0,
          nj=0,
          mj=0,
          Vj=0,
          Gj=0,
          Qj=0,
          Ij=0)
            if useSolution "To connect substance with solution, where is pressented"
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

         Chemical.Interfaces.SolutionState solutionState;
    equation

      if not useSolution then
        solutionState.T=solutionParam.T "Temperature of the solution";
        solutionState.p=solutionParam.p "Pressure of the solution";
        solutionState.v=solutionParam.v "Electric potential in the solution";
        solutionState.n=solutionParam.n "Amount of the solution";
        solutionState.m=solutionParam.m "Mass of the solution";
        solutionState.V=solutionParam.V "Volume of the solution";
        solutionState.G=solutionParam.G "Free Gibbs energy of the solution";
        solutionState.Q=solutionParam.Q "Electric charge of the solution";
        solutionState.I=solutionParam.I "Mole fraction based ionic strength of the solution";
      end if;

    end PartialSubstanceAndSolutionSource;
  end Internal;
annotation (Documentation(info="<html>
<u>This package contains the undirected junctions and necessary connectors between two undirected components, as well as directed-undirected connectors.</u>
<u>Note that in the undirected case it a distinction between junction and splitter is not possible.</u>
</html>", revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"), Icon(graphics={
        Line(
          points={{-80,2},{12,2}},
          color={158,66,200},
          thickness=0.5),
        Line(
          points={{12,2},{80,-78}},
          color={158,66,200},
          thickness=0.5),
        Line(
          points={{12,2},{80,82}},
          color={158,66,200},
          thickness=0.5),
        Ellipse(
          extent={{6,8},{18,-4}},
          lineColor={158,66,200},
          fillColor={194,138,221},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5)}));
end Topology;
