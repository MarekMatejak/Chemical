within Chemical;
package Topology "Junctions and Connectors for undirected chemical simulation"
  extends Modelica.Icons.Package;

  model JunctionMN "Generalized junction/splitter for undirected flow"

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Integer N(min=0) = 1 "Number of rears";
    parameter Integer M(min=0) = 1 "Number of fors";

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance for each branch" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of molar flow rate"
      annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rears[N](redeclare package stateOfMatter = stateOfMatter) "Rear ports"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld fores[M](redeclare package stateOfMatter = stateOfMatter) "Fore ports"
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
      fores[i].solution = rears[1].solution;

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

  model ConnectInletFore
    "Directed/undirected connector with input and fore"

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}),
          iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential r_fore;
    Modelica.Units.SI.ChemicalPotential r_inlet;

  equation
    fore.state_forwards = inlet.state;

    fore.n_flow + inlet.n_flow = 0;
    r_inlet + inlet.state.u =r_fore + Chemical.Utilities.Internal.regStep(
        inlet.n_flow,
        inlet.state.u,
        fore.state_rearwards.u,
        dropOfCommons.n_flow_reg);

    L/2*der(inlet.n_flow) = inlet.r - r_inlet;
    L/2*der(fore.n_flow) = fore.r - r_fore;

    fore.definition = inlet.definition;
    fore.solution = inlet.solution;

    annotation (Icon(
        graphics={
          Line(
            points={{-30,0},{30,0}},
            color={158,66,200},
            thickness=0.5)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect a unidirectional outlet to a undirected rear port. </u>
<u>The state from the inlet is given to the forward direction of the fore port, the total potential as well as the massflow of inlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
</html>"));
  end ConnectInletFore;

  model ConnectInletRear
    "Directed/undirected connector with input and rear"

   replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-40,0}),
          iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.RearOld rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={40,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
    ConnectInletFore connectInletFore(redeclare package stateOfMatter =
          stateOfMatter, final L=L/2)
      annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
    ConnectRearRear connectRearRear(redeclare package stateOfMatter =
          stateOfMatter, final L=L/2)
      annotation (Placement(transformation(extent={{2,-10},{22,10}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    connect(connectInletFore.fore, connectRearRear.rear_a) annotation (Line(
        points={{-7,0},{9,0}},
        color={158,66,200},
        thickness=0.5));
    connect(inlet, connectInletFore.inlet) annotation (Line(
        points={{-40,0},{-13,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectRearRear.rear_b, rear) annotation (Line(
        points={{15,0},{40,0}},
        color={158,66,200},
        thickness=0.5));
    annotation (Icon(
        graphics={
          Line(
            points={{-30,0},{30,0}},
            color={158,66,200},
            thickness=0.5)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect a unidirectional outlet to a undirected fore port. </u>
<u>The state from the inlet is given to the backward direction of the rear port, the total potential as well as the massflow of inlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
</html>"));
  end ConnectInletRear;

  model ConnectRearOutlet
    "Directed/undirected connector with rear and outlet"

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
    parameter Boolean useDefaultStateAsRear = false "Use Default Medium states as state_rearwards";

    Chemical.Interfaces.RearOld rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}),
          iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
    Chemical.Interfaces.StateInput state_rear if not useDefaultStateAsRear
      annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,40})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;
     Modelica.Units.SI.ChemicalPotential r_rear, r_outlet "Inertial potential";

  equation

    connect(rear.state_rearwards, state_rear);

    outlet.state = rear.state_forwards;

    rear.n_flow + outlet.n_flow = 0;
    r_outlet + outlet.state.u =r_rear + Chemical.Utilities.Internal.regStep(
        rear.n_flow,
        rear.state_forwards.u,
        rear.state_rearwards.u,
        dropOfCommons.n_flow_reg);

    L/2*der(outlet.n_flow) = outlet.r - r_outlet;
    L/2*der(rear.n_flow) = rear.r - r_rear;

    if useDefaultStateAsRear then
      rear.state_rearwards = Chemical.Interfaces.SubstanceState(0,0);
    end if;

    outlet.solution = rear.solution;
    outlet.definition = rear.definition;
    annotation (Icon(
        graphics={
          Line(
            points={{-30,0},{30,0}},
            color={158,66,200},
            thickness=0.5),
          Line(points={{2,58},{0,58}}, color={158,66,200}),
          Line(
            points={{0,0},{0,60}},
            color={162,29,33},
            arrow={Arrow.Filled,Arrow.None},
            arrowSize = 20)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect a unidirectional inlet to a undirected fore port. </u>
<u>The state from the forward direction of the rear port is handed to the outlet, the total potential as well as the massflow of outlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
<u>Since the unidirectional side of this connector is an oulet, there is no information about the rear state from this side. This information must be given by the StateInlet.</u>
</html>"));
  end ConnectRearOutlet;

  model ConnectForeOutlet
    "Directed/undirected connector with rear and outlet"
    extends
      Chemical.Topology.Internal.PartialSubstanceAndSolutionSource;

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
    parameter Boolean useDefaultStateAsRear = false "Use Default Medium states as state_rearwards";

    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-40,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={40,0}),
          iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));
    ConnectRearOutlet connectRearOutlet(redeclare package stateOfMatter =
          stateOfMatter, final L=L/2, final useDefaultStateAsRear = useDefaultStateAsRear)
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    ConnectForeFore connectForeFore(redeclare package stateOfMatter =
          stateOfMatter,
      substanceData=substanceData, solutionParam=solutionParam, useSolution=useSolution,
                         final L=L/2)
      annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
    Chemical.Interfaces.StateInput state_rear if not useDefaultStateAsRear
      annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,40})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation

    if useSolution then
      connect(solution,connectForeFore.solution);
    end if;
    connect(connectRearOutlet.rear, connectForeFore.fore_b) annotation (Line(
        points={{7,0},{-7,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectRearOutlet.outlet, outlet) annotation (Line(
        points={{13,0},{40,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectForeFore.fore_a, fore) annotation (Line(
        points={{-13,0},{-40,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectRearOutlet.state_rear, state_rear) annotation (Line(points={{10,4},{
            10,12},{0,12},{0,40}}, color={162,29,33}));
    annotation (Icon(
        graphics={
          Line(
            points={{-30,0},{30,0}},
            color={158,66,200},
            thickness=0.5), Line(
            points={{0,0},{0,60}},
            color={162,29,33},
            arrow={Arrow.Filled,Arrow.None},
            arrowSize = 20)},
        coordinateSystem(preserveAspectRatio=false)),
       Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>This connector can be used to connect a unidirectional inlet to a undirected rear port. </u>
<u>The state from the rearward direction of the fore port is handed to the outlet, the total potential as well as the massflow of outlet and port are set equal. </u>
<u>Note that when the flow is reversed, the resulting inertial potential can be different on both sides of this connector. </u>
<u>Since the unidirectional side of this connector is an oulet, there is no information about the rear state from this side. This information must be given by the StateInlet.</u>
</html>"));
  end ConnectForeOutlet;

  model ConnectForeFore "Undirected connector with fore and fore"
    extends Chemical.Topology.Internal.PartialSubstanceAndSolutionSource;
   /* replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby 
    Chemical.Interfaces.StateOfMatter
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
*/
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.ForeOld fore_a(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.ForeOld fore_b(redeclare package stateOfMatter = stateOfMatter)
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
    fore_b.solution=solutionState;
    fore_a.solution=solutionState;

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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rear_a(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={-30,0})));
    Chemical.Interfaces.RearOld rear_b(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={30,0}), iconTransformation(extent={{-20,-20},{20,20}}, origin={30,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    rear_b.state_rearwards = rear_a.state_forwards;
    rear_a.state_rearwards = rear_b.state_forwards;
    rear_a.n_flow + rear_b.n_flow = 0;
    L*der(rear_a.n_flow) = rear_a.r - rear_b.r;

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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld foreA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.ForeOld foreB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=2,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld foreA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.ForeOld foreB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=2,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
     parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.RearOld rearA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.RearOld rearB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=1,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
     parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    JunctionMN junctionMN(
      final M=1,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
      annotation (Placement(transformation(extent={{30,-10},{50,10}})));
    Chemical.Interfaces.RearOld rearA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    Chemical.Interfaces.RearOld rearB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld foreA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.ForeOld foreB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.ForeOld foreC(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=3,
      final N=1,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld reara(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.RearOld rearb(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.ForeOld foreB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    Chemical.Interfaces.ForeOld foreA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    JunctionMN junctionMN(
      final M=2,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld reara(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.RearOld rearb(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,-100})));
    Chemical.Interfaces.ForeOld foreA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.ForeOld foreB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=2,
      final N=2,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
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

    replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of each branch" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.RearOld rearA(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={0,100})));
    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,0})));
    JunctionMN junctionMN(
      final M=1,
      final N=3,
      final n_flow_reg=n_flow_reg,
      final L=L,
      redeclare package stateOfMatter=stateOfMatter)
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Chemical.Interfaces.RearOld rearB(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,0})));
    Chemical.Interfaces.RearOld rearC(redeclare package stateOfMatter = stateOfMatter)
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

  model ConnectorInletOutletFore

   replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
      Chemical.Interfaces.StateOfMatter
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

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_ref=dropOfCommons.n_flow_reg "Reference mass flow" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_ref=1e5 "Reference potential" annotation (Dialog(tab="Advanced"));
    parameter Boolean assumeConstantDensity = true "If true only mass-flow rate will determine the mixing"
      annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
      annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.ForeOld fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-20,-120},{20,-80}}), iconTransformation(extent={{-20,-120},{20,-80}})));
    Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{
              -80,20}})));
    Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
          stateOfMatter)
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,
              20}})));
    Chemical.FlowControl.CheckValve checkValve(
      redeclare package stateOfMatter = stateOfMatter,
      final L=L,
      final n_flow_ref=n_flow_ref,
      final u_ref=u_ref)
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
    Chemical.FlowControl.CheckValve checkValve1(
      redeclare package stateOfMatter = stateOfMatter,
      final L=L,
      final n_flow_ref=n_flow_ref,
      final u_ref=u_ref)
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));
    ConnectRearOutlet connectRearOutlet(
      redeclare package stateOfMatter = stateOfMatter,
      final L=L,
      final useDefaultStateAsRear=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    ConnectInletFore connectInletFore(
      redeclare package stateOfMatter = stateOfMatter,
      final L=L)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    JunctionRFF2 junctionRFF2_1(
      redeclare package stateOfMatter = stateOfMatter,
      final assumeConstantDensity=assumeConstantDensity,
      final n_flow_reg=n_flow_reg,
      final L=L)
        annotation (Placement(transformation(extent={{-10,10},{10,-10}})));
    Chemical.Sensors.SensorState sensorState(redeclare package stateOfMatter =
          stateOfMatter)
        annotation (Placement(transformation(extent={{-10,8},{10,28}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  equation
    connect(junctionRFF2_1.foreA, fore) annotation (Line(
        points={{0,-10},{0,-52},{0,-100},{0,-100}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF2_1.foreB, connectRearOutlet.rear)
      annotation (Line(
        points={{10,0},{27,0}},
        color={158,66,200},
        thickness=0.5));
    connect(junctionRFF2_1.rear, connectInletFore.fore) annotation (Line(
        points={{-10,0},{-27,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectInletFore.inlet, checkValve.outlet) annotation (Line(
        points={{-33,0},{-50,0}},
        color={158,66,200},
        thickness=0.5));
    connect(checkValve.inlet, inlet) annotation (Line(
        points={{-70,0},{-100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(connectRearOutlet.outlet, checkValve1.inlet) annotation (Line(
        points={{33,0},{50,0}},
        color={158,66,200},
        thickness=0.5));
    connect(checkValve1.outlet, outlet) annotation (Line(
        points={{70,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(sensorState.inlet, checkValve.outlet)
      annotation (Line(
        points={{-10,18},{-40,18},{-40,0},{-50,0}},
        color={158,66,200},
        thickness=0.5));
    connect(sensorState.state_out, connectRearOutlet.state_rear) annotation (Line(points={{10,18},{30,18},{30,4}}, color={162,29,33}));
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
            points={{0,-100},{0,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,6},{6,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConnectorInletOutletFore;

  package Tests "Test models for the undirected topology package"
  extends Modelica.Icons.ExamplesPackage;

    model TestJunction "Test for the undirected junction"
      extends Modelica.Icons.Example;

      replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medum model for the Test" annotation (Documentation(info="<html>
<u>This is the replaceable package that determines the medium of the Test. </u>
</html>"));

      replaceable function pLoss =
          Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss
          (
        material=Chemical.Processes.Internal.Material.wood)
        constrainedby Chemical.Processes.Internal.FlowResistance.partialChemicalPotentialLoss
        "ChemicalPotential loss function for all Flow resistances"
        annotation(choicesAllMatching = true, Documentation(info="<html>
<u>This is the potential loss function used for all resistances except the two on the outlets of the right two cases.</u>
</html>"));

      Chemical.Boundaries.Source source(redeclare package Medium=Medium,
          u0_par=200000)
        annotation (Placement(transformation(extent={{-160,30},{-140,50}})));
      Chemical.Boundaries.Sink sink(redeclare package Medium=Medium,
          u0_par=150000)
        annotation (Placement(transformation(extent={{-20,50},{0,70}})));
      Chemical.Boundaries.Sink sink1(redeclare package Medium=Medium,
          u0_par=100000)
        annotation (Placement(transformation(extent={{-20,10},{0,30}})));
      Chemical.Topology.SplitterT1 splitterT1_1(redeclare package Medium=Medium, L=
            2*dropOfCommons.L)
        annotation (Placement(transformation(extent={{-94,30},{-74,50}})));
      inner Chemical.DropOfCommons dropOfCommons(L=1e2) annotation (Placement(transformation(extent={{-134,-8},{-114,12}})));
      Chemical.Processes.FlowResistance flowResistance(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-46,50},{-26,70}})));
      Chemical.Processes.FlowResistance flowResistance1(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-46,10},{-26,30}})));
      Chemical.Boundaries.Source source2(
        redeclare package stateOfMatter = stateOfMatter, u0_par=2000000)
        annotation (Placement(transformation(extent={{0,50},{20,70}})));
      Chemical.Boundaries.Source source3(
        redeclare package stateOfMatter = stateOfMatter, u0_par=3500000)
        annotation (Placement(transformation(extent={{0,10},{20,30}})));
      Chemical.Processes.FlowResistance flowResistance4(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{26,10},{46,30}})));
      Chemical.Processes.FlowResistance flowResistance5(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{26,50},{46,70}})));
      Chemical.Topology.JunctionT1 junctionT1_1(redeclare package Medium=Medium,
        assumeConstantDensity=false,
        L=2*dropOfCommons.L)
        annotation (Placement(transformation(extent={{94,30},{74,50}})));
      Chemical.Boundaries.Sink sink4(
        redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{140,30},{160,50}})));
      Chemical.Boundaries.Source source6(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=200000)
        annotation (Placement(transformation(extent={{-160,-50},{-140,-30}})));
      Chemical.Boundaries.Sink sink6(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=150000)
        annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
      Chemical.Boundaries.Sink sink7(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=100000)
        annotation (Placement(transformation(extent={{-20,-70},{0,-50}})));
      ConnectInletFore connectInletFore1(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-118,-50},{-98,-30}})));
      ConnectRearOutlet connectRearOutlet2(redeclare package stateOfMatter =
            stateOfMatter,
          useDefaultStateAsRear=true)
        annotation (Placement(transformation(extent={{-66,-30},{-46,-10}})));
      ConnectRearOutlet connectRearOutlet3(redeclare package stateOfMatter =
            stateOfMatter,
          useDefaultStateAsRear=true)
        annotation (Placement(transformation(extent={{-66,-70},{-46,-50}})));
      Chemical.Processes.FlowResistance flowResistance8(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-46,-30},{-26,-10}})));
      Chemical.Processes.FlowResistance flowResistance9(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-46,-70},{-26,-50}})));
      Chemical.Boundaries.Source source7(redeclare package stateOfMatter =
            stateOfMatter,                                                                u0_par=2000000)
        annotation (Placement(transformation(extent={{0,-30},{20,-10}})));
      Chemical.Boundaries.Source source8(redeclare package stateOfMatter =
            stateOfMatter,                                                                u0_par=3500000)
        annotation (Placement(transformation(extent={{0,-70},{20,-50}})));
      Chemical.Processes.FlowResistance flowResistance10(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{26,-70},{46,-50}})));
      Chemical.Processes.FlowResistance flowResistance11(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{26,-30},{46,-10}})));
      Chemical.Boundaries.Sink sink8(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=100000)
        annotation (Placement(transformation(extent={{140,-50},{160,-30}})));
      ConnectInletRear connectInletRear2(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{46,-30},{66,-10}})));
      ConnectInletRear connectInletRear3(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{46,-70},{66,-50}})));
      ConnectForeOutlet connectForeOutlet1(redeclare package stateOfMatter =
            stateOfMatter,
          useDefaultStateAsRear=true)
        annotation (Placement(transformation(extent={{98,-50},{118,-30}})));
      JunctionMN junctionMN(redeclare package stateOfMatter = stateOfMatter,
        M=2,
        N=1)
        annotation (Placement(transformation(extent={{-94,-50},{-74,-30}})));
      JunctionMN junctionMN1(redeclare package stateOfMatter = stateOfMatter,
        M=2,
        N=1,
        assumeConstantDensity=false)
        annotation (Placement(transformation(extent={{94,-50},{74,-30}})));
      Chemical.Processes.FlowResistance flowResistance2(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-136,30},{-116,50}})));
      Chemical.Processes.FlowResistance flowResistance3(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = pLoss)
        annotation (Placement(transformation(extent={{-136,-50},{-116,-30}})));
      Chemical.Processes.FlowResistance flowResistance6(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
            (k=1000))
        annotation (Placement(transformation(extent={{116,-50},{136,-30}})));
      Chemical.Processes.FlowResistance flowResistance7(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
            (k=1000))
        annotation (Placement(transformation(extent={{116,30},{136,50}})));
    equation
      connect(sink.inlet, flowResistance.outlet) annotation (Line(
          points={{-20,60},{-26,60}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance.inlet, splitterT1_1.outletA) annotation (Line(
          points={{-46,60},{-84,60},{-84,50}},
          color={158,66,200},
          thickness=0.5));
      connect(sink1.inlet, flowResistance1.outlet) annotation (Line(
          points={{-20,20},{-26,20}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance1.inlet, splitterT1_1.outletB) annotation (Line(
          points={{-46,20},{-84,20},{-84,30}},
          color={158,66,200},
          thickness=0.5));
      connect(source2.outlet, flowResistance5.inlet) annotation (Line(
          points={{20,60},{26,60}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance4.inlet, source3.outlet) annotation (Line(
          points={{26,20},{20,20}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance4.outlet, junctionT1_1.inletB) annotation (Line(
          points={{46,20},{84,20},{84,30}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionT1_1.inletA, flowResistance5.outlet) annotation (Line(
          points={{84,50},{84,60},{46,60}},
          color={158,66,200},
          thickness=0.5));
      connect(sink6.inlet,flowResistance8. outlet) annotation (Line(
          points={{-20,-20},{-26,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance8.inlet, connectRearOutlet2.outlet) annotation (Line(
          points={{-46,-20},{-53,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance9.inlet,connectRearOutlet3. outlet) annotation (Line(
          points={{-46,-60},{-53,-60}},
          color={158,66,200},
          thickness=0.5));
      connect(sink7.inlet,flowResistance9. outlet) annotation (Line(
          points={{-20,-60},{-26,-60}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance10.inlet, source8.outlet) annotation (Line(
          points={{26,-60},{20,-60}},
          color={158,66,200},
          thickness=0.5));
      connect(source7.outlet, flowResistance11.inlet) annotation (Line(
          points={{20,-20},{26,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletRear2.inlet, flowResistance11.outlet) annotation (Line(
          points={{53,-20},{46,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletRear3.inlet, flowResistance10.outlet) annotation (Line(
          points={{53,-60},{46,-60}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet2.rear, junctionMN.fores[2]) annotation (Line(
          points={{-59,-20},{-66,-20},{-66,-39.5},{-74,-39.5}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet3.rear, junctionMN.fores[1]) annotation (Line(
          points={{-59,-60},{-66,-60},{-66,-40.5},{-74,-40.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN.rears[1], connectInletFore1.fore) annotation (Line(
          points={{-94,-40},{-100,-40},{-100,-40},{-105,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionMN1.rears[1], connectForeOutlet1.fore) annotation (Line(
          points={{94,-40},{100,-40},{100,-40},{105,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletRear2.rear, junctionMN1.fores[2]) annotation (Line(
          points={{59,-20},{66,-20},{66,-39.5},{74,-39.5}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletRear3.rear, junctionMN1.fores[1]) annotation (Line(
          points={{59,-60},{66,-60},{66,-40.5},{74,-40.5}},
          color={158,66,200},
          thickness=0.5));
      connect(source.outlet, flowResistance2.inlet) annotation (Line(
          points={{-140,40},{-136,40}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterT1_1.inlet, flowResistance2.outlet) annotation (Line(
          points={{-94,40},{-116,40}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletFore1.inlet, flowResistance3.outlet) annotation (Line(
          points={{-111,-40},{-116,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(source6.outlet, flowResistance3.inlet) annotation (Line(
          points={{-140,-40},{-136,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(connectForeOutlet1.outlet, flowResistance6.inlet) annotation (Line(
          points={{111,-40},{116,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(sink8.inlet, flowResistance6.outlet) annotation (Line(
          points={{140,-40},{136,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(sink4.inlet, flowResistance7.outlet) annotation (Line(
          points={{140,40},{136,40}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionT1_1.outlet, flowResistance7.inlet) annotation (Line(
          points={{94,40},{116,40}},
          color={158,66,200},
          thickness=0.5));
      annotation (
      experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
      Documentation(info="<html>
<u>This model tests the undirected junction against the unidirectional junction. The states of the left two and the right two systems are expected to be the same, when the junctions have the same settings. </u>
<u>Note that the unidirectional junctions have two times the inertance, since the undirected junction comes with a additional connector, which in turn adds inertance to each outlet.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"),
        Diagram(coordinateSystem(extent={{-160,-80},{160,80}})),
        Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
    end TestJunction;

    model TestConnectors "Test for the connectors"
      extends Modelica.Icons.Example;

      replaceable package Medium = Chemical.Media.myMedia.Air.SimpleAir constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
                                                             "Medum model for the Test"
        annotation (Documentation(info="<html>
<u>This is the replaceable package that determines the medium of the Test. </u>
</html>"));

      ConnectRearRear connectRearRear(redeclare package Medium=Medium)
        annotation (Placement(transformation(extent={{-10,60},{10,80}})));
      ConnectForeFore connectForeFore(redeclare package Medium=Medium)
        annotation (Placement(transformation(extent={{-10,40},{10,60}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
        redeclare package stateOfMatter = stateOfMatter,
        u0_par=100000,
        fore(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
      Chemical.Boundaries.BoundaryFore boundary_fore(redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{-20,40},{-40,60}})));
      Chemical.Boundaries.Source source(redeclare package Medium=Medium,
          u0_par=100000)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
      ConnectInletFore connectInletFore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      ConnectRearOutlet connectRearOutlet(redeclare package Medium=Medium, rear(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-10,0},{10,20}})));
      Chemical.Boundaries.Sink sink(redeclare package Medium=Medium,
          potentialFromInput=true)
        annotation (Placement(transformation(extent={{20,0},{40,20}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(redeclare package stateOfMatter = stateOfMatter, potentialFromInput=true)
        annotation (Placement(transformation(extent={{40,60},{20,80}})));
      Chemical.Boundaries.BoundaryFore boundary_fore1(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,40},{40,60}})));
      Chemical.Boundaries.BoundaryFore boundary_fore2(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,20},{40,40}})));
      Chemical.Boundaries.BoundaryRear boundary_rear2(redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
      Chemical.Boundaries.Source source1(redeclare package stateOfMatter =
            stateOfMatter,
        u0_par=100000,
        L=1.5*dropOfCommons.L,
        outlet(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-40,100},{-20,120}})));
      Chemical.Boundaries.Sink sink1(redeclare package stateOfMatter =
            stateOfMatter,
        potentialFromInput=true,
        L=1.5*dropOfCommons.L)
        annotation (Placement(transformation(extent={{20,100},{40,120}})));
      Chemical.Boundaries.BoundaryRear boundary_rear3(
        redeclare package stateOfMatter = stateOfMatter,
        u0_par=100000,
        L=1.5*dropOfCommons.L) annotation (Placement(transformation(extent={{-40,80},{-20,100}})));
      Chemical.Boundaries.BoundaryFore boundary_fore3(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        L=1.5*dropOfCommons.L,
        rear(n_flow(start=0, fixed=true))) annotation (Placement(transformation(extent={{20,80},{40,100}})));
      inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{-88,80},{-68,100}})));
      Modelica.Blocks.Sources.Pulse pulse(
        amplitude=1e5,
        period=5,
        offset=0.5e5,
        startTime=1)
        annotation (Placement(transformation(extent={{82,50},{62,70}})));
      Chemical.Boundaries.Source source2(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=100000)
        annotation (Placement(transformation(extent={{-50,-40},{-30,-20}})));
      ConnectInletFore connectInletFore1(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{-30,-40},{-10,-20}})));
      ConnectRearOutlet connectRearOutlet1(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(extent={{10,-40},{30,-20}})));
      Chemical.Boundaries.Sink sink2(redeclare package stateOfMatter =
            stateOfMatter,
          potentialFromInput=true)
        annotation (Placement(transformation(extent={{30,-40},{50,-20}})));
      Chemical.Processes.FlowResistance flowResistance(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        computeL=false,
        r=0.01,
        l=1,
        redeclare function pLoss = .Chemical.Processes.Internal.FlowResistance.laminarChemicalPotentialLoss)
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
      Chemical.Processes.FlowResistance flowResistance2(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        computeL=false,
        r=0.01,
        l=1,
        redeclare function pLoss = .Chemical.Processes.Internal.FlowResistance.laminarChemicalPotentialLoss)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
      Chemical.Boundaries.BoundaryRear boundary_rear4(redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
      Chemical.Boundaries.BoundaryFore boundary_fore4(redeclare package stateOfMatter = stateOfMatter, potentialFromInput=true)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={40,-50})));
      Chemical.Boundaries.CreateState createState(
        redeclare package stateOfMatter = stateOfMatter,
        PFromInput=true)
        annotation (Placement(transformation(extent={{38,-12},{32,-4}})));
      ConnectorInletOutletFore connectorInletOutletFore(redeclare package stateOfMatter =
            stateOfMatter)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,-68})));
      Chemical.Boundaries.BoundaryFore boundary_fore5(redeclare package stateOfMatter = stateOfMatter, potentialFromInput=true)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={64,-90})));
      Chemical.Boundaries.Source source3(redeclare package stateOfMatter =
            stateOfMatter,                                                                u0_par=200000)
        annotation (Placement(transformation(extent={{-80,-78},{-60,-58}})));
      Chemical.Boundaries.Sink sink3(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=false,
        u0_par=100000)
        annotation (Placement(transformation(extent={{54,-78},{74,-58}})));
      Chemical.Processes.FlowResistance flowResistance1(
        redeclare package stateOfMatter = stateOfMatter,
        r(displayUnit="mm") = 0.005,
        l=5,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss)
        annotation (Placement(transformation(extent={{-40,-78},{-20,-58}})));
      Chemical.Processes.FlowResistance flowResistance3(
        redeclare package stateOfMatter = stateOfMatter, initM_flow = Chemical.Utilities.Types.InitializationMethods.state,
        l=5,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss,
        r(displayUnit="mm") = 0.005)
        annotation (Placement(transformation(extent={{20,-78},{40,-58}})));
      Chemical.Processes.FlowResistance flowResistance4(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        l=5,
        redeclare function pLoss = Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss,
        r(displayUnit="mm") = 0.005) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={30,-90})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=1.4e5,
        duration=1,
        offset=0.8e5,
        startTime=4.5) annotation (Placement(transformation(extent={{100,-100},{80,-80}})));
    equation
      connect(boundary_rear.fore, connectRearRear.rear_a) annotation (Line(
          points={{-20,70},{-3,70}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore.rear, connectForeFore.fore_a) annotation (Line(
          points={{-20,50},{-3,50}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletFore.inlet, source.outlet) annotation (Line(
          points={{-3,30},{-20,30}},
          color={158,66,200},
          thickness=0.5));
      connect(sink.inlet, connectRearOutlet.outlet) annotation (Line(
          points={{20,10},{3,10}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear1.fore, connectRearRear.rear_b) annotation (Line(
          points={{20,70},{3,70}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore1.rear, connectForeFore.fore_b) annotation (Line(
          points={{20,50},{3,50}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore2.rear, connectInletFore.fore) annotation (Line(
          points={{20,30},{3,30}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear2.fore, connectRearOutlet.rear) annotation (Line(
          points={{-20,10},{-3,10}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear3.fore, boundary_fore3.rear) annotation (Line(
          points={{-20,90},{20,90}},
          color={158,66,200},
          thickness=0.5));
      connect(sink1.inlet, source1.outlet) annotation (Line(
          points={{20,110},{-20,110}},
          color={158,66,200},
          thickness=0.5));
      connect(pulse.y, sink1.u0_var) annotation (Line(points={{61,60},{46,60},{46,110},{32,110}},
                                                                                              color={0,0,127}));
      connect(pulse.y, boundary_fore3.u0_var) annotation (Line(points={{61,60},{46,60},{46,96},{32,96}},
                                                                                                       color={0,0,127}));
      connect(pulse.y, boundary_rear1.u0_var) annotation (Line(points={{61,60},{46,60},{46,76},{32,76}},
                                                                                                     color={0,0,127}));
      connect(pulse.y, boundary_fore1.u0_var) annotation (Line(points={{61,60},{46,60},{46,56},{32,56}}, color={0,0,127}));
      connect(pulse.y, boundary_fore2.u0_var) annotation (Line(points={{61,60},{46,60},{46,36},{32,36}}, color={0,0,127}));
      connect(pulse.y, sink.u0_var) annotation (Line(points={{61,60},{46,60},{46,10},{32,10}}, color={0,0,127}));
      connect(connectInletFore1.inlet, source2.outlet) annotation (Line(
          points={{-23,-30},{-30,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(sink2.inlet, connectRearOutlet1.outlet) annotation (Line(
          points={{30,-30},{23,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(connectInletFore1.fore, flowResistance.rear) annotation (Line(
          points={{-17,-30},{-10,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(connectRearOutlet1.rear, flowResistance.fore) annotation (Line(
          points={{17,-30},{10,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(sink2.u0_var, sink.u0_var) annotation (Line(points={{42,-30},{46,-30},{46,10},{32,10}},
                                  color={0,0,127}));
      connect(flowResistance2.rear, boundary_rear4.fore) annotation (Line(
          points={{-10,-50},{-30,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance2.fore, boundary_fore4.rear) annotation (Line(
          points={{10,-50},{30,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore4.u0_var, sink.u0_var) annotation (Line(points={{42,-44},{46,-44},{46,10},{32,10}},
                                            color={0,0,127}));
      connect(createState.y, connectRearOutlet.state_rear) annotation (Line(points={{32,-8},{12,-8},{12,20},{0,20},{0,14}},
                                                                            color={
              162,29,33}));
      connect(createState.u_inp, sink.u0_var) annotation (Line(points={{38,-4},{46,-4},{46,10},{32,10}},
                                       color={0,0,127}));
      connect(connectRearOutlet1.state_rear, connectRearOutlet.state_rear)
        annotation (Line(points={{20,-26},{20,-8},{12,-8},{12,28},{0,28},{0,14}},
            color={162,29,33}));
      connect(source3.outlet, flowResistance1.inlet) annotation (Line(
          points={{-60,-68},{-40,-68}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance1.outlet, connectorInletOutletFore.inlet)
        annotation (Line(
          points={{-20,-68},{-10,-68}},
          color={158,66,200},
          thickness=0.5));
      connect(sink3.inlet, flowResistance3.outlet) annotation (Line(
          points={{54,-68},{40,-68}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance3.inlet, connectorInletOutletFore.outlet)
        annotation (Line(
          points={{20,-68},{10,-68}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore5.rear, flowResistance4.fore) annotation (Line(
          points={{54,-90},{40,-90}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance4.rear, connectorInletOutletFore.fore)
        annotation (Line(
          points={{20,-90},{0,-90},{0,-78}},
          color={158,66,200},
          thickness=0.5));
      connect(ramp.y, boundary_fore5.u0_var) annotation (Line(points={{79,-90},{72,-90},{72,-84},{66,-84}},
                                                                                          color={0,0,127}));
      annotation (
      experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
      Documentation(info="<html>
<u>This model tests the connectors against the case when source is directly connected to sink. All massflows are expected to be the same, when the sources and sinks have the same settings. </u>
<u>Note that the sources and sinks without a connector have 1.5 times the inertance, since the additional connector adds inertance to each other path.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"),
        Diagram(coordinateSystem(extent={{-100,-140},{100,120}})),
        Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
    end TestConnectors;
    annotation (Documentation(info="<html>
<u>This package contains the test cases for the undirected topology package. </u>
</html>"));
  end Tests;

  package Internal
    partial model PartialSubstanceAndSolutionSource
      "Substance properties for components, where the substance is connected with the solution"

     replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby
        Chemical.Interfaces.StateOfMatter
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

      parameter stateOfMatter.SubstanceDataParameters substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable=not useRear));

      parameter Boolean useSolution = false "Use solution connector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
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
