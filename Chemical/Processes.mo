within Chemical;
package Processes "Undirected process package"
  model Reaction "Chemical Reaction"
    extends Chemical.Processes.Internal.PartialReactionWithSubstanceData;
    extends Chemical.Interfaces.ConditionalKinetics
                                          (k_forward=1);
    Real rr_fore_exact,rr_rear_exact,  kb;
  equation

    rr = kf * Sx_fore * ( 1  -  exp(-duRT_fore));
    if nS>0 then
      rr = -kb * Px_rear * ( 1  -  exp(duRT_rear));
    end if;

    Kx = kb/kf;

    //the same as:
    rr_fore_exact = (kf*Sx_fore - kb*Px_fore);
    rr_rear_exact = (kf*Sx_rear - kb*Px_rear);

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Rectangle(
            extent={{-100,-30},{100,30}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-100,-72},{100,-40}},
            lineColor={128,0,255},
          textString="%name"),
          Polygon(
            points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
<td><p>molar Gibb&apos;s energy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>H = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>H<sub>i</sub>) </p></td>
<td><p>molar enthalpy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>S = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>S<sub>i</sub>) = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>r</sub>&omega;) </p></td>
<td><p>molar entropy of the reaction</p></td>
</tr>
</table>
<h4><span style=\"color:#008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stochiometric coefficients of i-th substance</p></td>
</tr>
<tr>
<td><p>K</p></td>
<td><p>dissociation constant (activity based)</p></td>
</tr>
<tr>
<td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
<td><p>activity of the substance A</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub></p></td>
<td><p>activity coefficient of the substance A</p></td>
</tr>
<tr>
<td><p>x<sub>i</sub></p></td>
<td><p>mole fraction of the substance A</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
<td><p>molar enthalpy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
<td><p>molar Gibbs energy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
<td><p>molar entropy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by reaction</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
  end Reaction;
  extends Modelica.Icons.Package;

  model ForwardReaction "Chemical Reaction"
    extends Internal.PartialReactionWithSubstanceData;
    extends Chemical.Interfaces.ConditionalKinetics
                                          (k_forward=1);

  equation

    rr = kf * Sx_fore;

    if nP>0 then
       der(products[1].state_forwards.u).*TC = products[1].r;
    end if;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Rectangle(
            extent={{-100,-30},{100,30}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-100,-72},{100,-40}},
            lineColor={128,0,255},
          textString="%name"),
          Polygon(
            points={{-60,2},{-60,0},{54,0},{54,0},{18,10},{18,2},{-60,2}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,-2},{-60,0},{54,0},{54,0},{18,-10},{18,-2},{-60,-2}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
<td><p>molar Gibb&apos;s energy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>H = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>H<sub>i</sub>) </p></td>
<td><p>molar enthalpy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>S = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>S<sub>i</sub>) = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>r</sub>&omega;) </p></td>
<td><p>molar entropy of the reaction</p></td>
</tr>
</table>
<h4><span style=\"color:#008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stochiometric coefficients of i-th substance</p></td>
</tr>
<tr>
<td><p>K</p></td>
<td><p>dissociation constant (activity based)</p></td>
</tr>
<tr>
<td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
<td><p>activity of the substance A</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub></p></td>
<td><p>activity coefficient of the substance A</p></td>
</tr>
<tr>
<td><p>x<sub>i</sub></p></td>
<td><p>mole fraction of the substance A</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
<td><p>molar enthalpy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
<td><p>molar Gibbs energy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
<td><p>molar entropy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by reaction</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
  end ForwardReaction;

  model FlowResistance "Flow resistance model"
    extends Chemical.Interfaces.SISO                 (
                                  final L=if computeL then l/(r^2*pi) else L_value, final cliu_u_out=true);

    import Modelica.Constants.pi "Constant Pi";

    parameter Modelica.Units.SI.Radius r(min=0) "Radius of pipe";
    parameter Modelica.Units.SI.Length l(min=0) "Length of pipe";
    parameter Chemical.Utilities.Units.Inertance L_value=dropOfCommons.L "Inertance of pipe" annotation (Dialog(enable=not computeL, tab="Advanced"));
    parameter Boolean computeL = true "Compute L from r and l"
      annotation(Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal input density" annotation (Dialog(tab="Advanced"));

    replaceable function pLoss =
        Chemical.Processes.Internal.FlowResistance.pleaseSelectChemicalPotentialLoss
      constrainedby
      Chemical.Processes.Internal.FlowResistance.partialChemicalPotentialLoss "ChemicalPotential loss function"
      annotation(choicesAllMatching=true, Documentation(info="<html>
<u>ChemicalPotential loss function used in the flow resistance.</u>
</html>"));

  protected
    Modelica.Units.SI.Density rho_rear_in=max(rho_min, Medium.density(rear.state_forwards)) "density of medium entering";
    Modelica.Units.SI.DynamicViscosity mu_rear_in=Medium.dynamicViscosity(rear.state_forwards) "dynamic viscosity of medium entering";

    Modelica.Units.SI.Density rho_fore_in=max(rho_min, Medium.density(fore.state_rearwards)) "density of medium entering";
    Modelica.Units.SI.DynamicViscosity mu_fore_in=Medium.dynamicViscosity(fore.state_rearwards) "dynamic viscosity of medium entering";

  equation
    //forwards model
    du_fore = -pLoss(n_flow, rho_rear_in, mu_rear_in, r, l);
    h_fore_out = h_rear_in;
    Xi_fore_out = Xi_rear_in;

    //rearwards model
    du_rear = -pLoss(-n_flow, rho_fore_in, mu_fore_in, r, l);
    h_rear_out = h_fore_in;
    Xi_rear_out = Xi_fore_in;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Ellipse(
            extent={{-56,54},{64,-66}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-100,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-60,60},{60,-60}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{40,0},{-48,0}},
            color={158,66,200},
            thickness=0.5,
            pattern=LinePattern.Dash),
          Line(
            points={{-44,-40},{0,-10},{44,-40}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{-44,-15},{0,15},{44,-15}},
            color={158,66,200},
            thickness=0.5,
            smooth=Smooth.Bezier,
            origin={0,25},
            rotation=180)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Undirected implementation of the FlowResistance with different selectable flow resistance functions (laminar, laminar-turbulent, linear-quadratic). The output potential can be clipped to a certain value.</u>
</html>"));
  end FlowResistance;

  model TransportDelay "Delay chemical state depending on fluid speed"
    extends Chemical.Interfaces.SISO                 (final cliu_u_out=
          false);

    parameter Modelica.Units.SI.Length l "Length of Delay Pipe";
    parameter Modelica.Units.SI.Radius r "Radius of Delay Pipe";
    parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal Density" annotation (Dialog(tab="Advanced"));
    parameter Real v_min(unit="1/s", min=0) = 0.01 "Minimum nondimensional speed"
      annotation(Dialog(tab="Advanced"));
    parameter Real v_max(unit="1/s", min=0) = 50 "Maximum nondimensional speed"
      annotation(Dialog(tab="Advanced"));

    constant Medium.ThermodynamicState state_0 = Chemical.Interfaces.SubstanceState(u=Medium.u_default,h= Medium.h_default);

    constant Modelica.Units.SI.SpecificVolume v_0=1/Medium.density(state_0);

    Real x(unit="1");
    Real v(unit="1/s");

  protected
    Modelica.Units.SI.SpecificInternalEnergy u_rear_in=Medium.specificInternalEnergy(rear.state_forwards);
    Modelica.Units.SI.SpecificInternalEnergy u_fore_in=Medium.specificInternalEnergy(fore.state_rearwards);
    Modelica.Units.SI.SpecificVolume v_rear_in=1/max(rho_min, Medium.density(rear.state_forwards));
    Modelica.Units.SI.SpecificVolume v_fore_in=1/max(rho_min, Medium.density(fore.state_rearwards));

    Modelica.Units.SI.SpecificInternalEnergy u_rear_out;
    Modelica.Units.SI.SpecificInternalEnergy u_fore_out;
    Modelica.Units.SI.SpecificVolume v_rear_out;
    Modelica.Units.SI.SpecificVolume v_fore_out;

    Modelica.Units.SI.Area A=r^2*Modelica.Constants.pi;

  initial equation
    x = 0;

  equation
    if n_flow >= 0 then
      v = min(v_max, max(v_min, n_flow*v_rear_in/A/l));
    else
      v = -min(v_max, max(v_min, -n_flow*v_fore_in/A/l));
    end if;

    der(x) = v;

    (u_rear_out,u_fore_out) = spatialDistribution(u_rear_in, u_fore_in,
      x, v>=0,
      initialPoints = {0.0,1.0},
      initialValues = {Medium.specificInternalEnergy(state_0), Medium.specificInternalEnergy(state_0)});
    (v_rear_out,v_fore_out) = spatialDistribution(v_rear_in, v_fore_in,
      x, v>=0,
      initialPoints = {0.0,1.0},
      initialValues = {v_0, v_0});

    for i in 1:Medium.nXi loop
      (Xi_rear_out[i], Xi_fore_out[i]) = spatialDistribution(Xi_rear_in[i], Xi_fore_in[i],
        x, v>=0,
        initialPoints = {0.0,0.0},
        initialValues = {Xi_0[i], Xi_0[i]});
    end for;

      //forwards model
    du_fore = 0;
    h_fore_out = u_fore_out + u_fore_out * v_fore_out;
    Xi_fore_out = Xi_rear_in;

    //rearwards model
    du_rear = 0;
    h_rear_out = u_rear_out + u_rear_out * v_rear_out;
    Xi_rear_out = Xi_fore_in;
    annotation (Documentation(info="<html>
<u>Undirected implementation of the transport delay.</u>
<u>Delays the temperature and massFraction, not potential, since potential differences propagate with speed of sound, and since delaying only steady state potential u not inertial potential r might lead to undesirable behavior.</u>
<u>Note that this component uses the spatialDistribution operator, that has some artefacts (see Fig. 1) for high and low non-dimensional speeds v (possibly due to inerpolation or extrapolation of the function). Therefore minimum and maximum speed in the non-dimensional coordinate x (inlet @ x=0, outlet @ x=1) is limited. The default limits are [0.01, 50], so the delay is limited by default to [0.02s, 100s]. This limit can be adjusted in the advanced parameters tab.</u>
<u><img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.Tests.TransportDelay_artefacts2.PNG\"/> <img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.Tests.TransportDelay_artefacts.PNG\"/> </u>
<u style=\"margin-left: 250px;\">Fig. 1: artefacts of the TransportDelay</u>
</html>"),
  Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Ellipse(
            extent={{-56,54},{64,-66}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-100,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-60,60},{60,-60}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
             points={{-25,-43},{25,-43},{-25,43},{25,43}},
             lineColor= {158,66,200},
             lineThickness=0.5),
          Line(
            points={{0,0},{0,-30}},
            color={158,66,200},
            thickness=0.5,
            pattern=LinePattern.Dot),
          Polygon(
             points={{-25,-43},{25,-43}, {0,-30}},
             lineColor= {158,66,200},
             fillColor={158,66,200},
             fillPattern=FillPattern.Solid),
          Polygon(
             points={{-15,26},{15,26}, {0,0}},
             lineColor= {158,66,200},
             fillColor={158,66,200},
             fillPattern=FillPattern.Solid)}));
  end TransportDelay;

  model ConductionElement "Volume with quasi-sationary mass and heatport"
    extends Internal.PartialConductionElement;

    parameter Boolean resistanceFromAU = true
      "= true, if thermal conductance given by U*A"
      annotation(Dialog(group="Thermal Conductance"));
    parameter Modelica.Units.SI.Area A=1 "Contact area of element with medium" annotation (Dialog(group="Thermal Conductance", enable=resistanceFromAU));
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer U=200 "Heat transfer coefficient to medium"
      annotation (Dialog(group="Thermal Conductance", enable=resistanceFromAU));
    parameter Modelica.Units.SI.ThermalConductance k_par=200 "Thermal conductance heatport->fluid"
      annotation (Dialog(group="Thermal Conductance", enable=not resistanceFromAU));

  equation
    k = (if noEvent(resistanceFromAU) then U*A else k_par);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>
Undirected implementation of the Conduction Element.
</u>
<u>
This model is an element with a fixed volume (fig. 1). The mass in the volume is
assumed quasi-stationary (statically computed with volume and density), and the
fore massflow is coupled to the rear massflow.
</u>
<u>
<strong>Because of this the ConductionElement cannot be used as a loop breaker.</strong>
</u>
<u>
The advantage is that multiple ConductionElements can be put behind each other
without worrying about oscillations or fast eigenvalues between their masses.
The ConductionElement implements equations for conservation of mass and energy
for the fluid mass contained within it.
</u>
<u>
For further documentation see the documentation of the
<a href=\"modelica://Chemical.Undirected.Processes.Internal.PartialConductionElement\">motherclass</a>.
</u>
</html>"));
  end ConductionElement;

  package Tests "Tests for top level components of undirected"
    extends Modelica.Icons.ExamplesPackage;

    model TestFlowResistance "Test for the undirected flow resistance"
      extends Modelica.Icons.Example;

      Chemical.Boundaries.BoundaryRear boundary_rear(substanceData=Chemical.Substances.Water_liquid(), u0_par=100000)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-30,0})));
      Chemical.Boundaries.BoundaryFore boundary_fore(potentialFromInput=true, u0_par=110000) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      Modelica.Blocks.Sources.Step step(
        height=-80000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{60,-6},{48,6}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(substanceData=Chemical.Substances.Water_liquid(), potentialFromInput=true)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-28,-38})));
      Chemical.Boundaries.BoundaryFore boundary_fore1(potentialFromInput=false, u0_par=100000)
        annotation (Placement(transformation(extent={{22,-48},{42,-28}})));
      Reaction reaction(
        productsSubstanceData={Chemical.Substances.Water_liquid()},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Reaction reaction1(
        productsSubstanceData={Chemical.Substances.Water_liquid()},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,-48},{12,-28}})));
    equation
      connect(step.y, boundary_fore.u0_var)
        annotation (Line(points={{47.4,0},{40,0},{40,6},{32,6}},
                                                       color={0,0,127}));
      connect(boundary_rear1.u0_var, boundary_fore.u0_var)
        annotation (Line(points={{-30,-44},{-40,-44},{-40,-20},{40,-20},{40,6},{32,6}}, color={0,0,127}));
      connect(boundary_rear.fore, reaction.substrates[1]) annotation (Line(
          points={{-20,0},{-10,0}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction.products[1], boundary_fore.rear) annotation (Line(
          points={{10,0},{20,0}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear1.fore, reaction1.substrates[1]) annotation (Line(
          points={{-18,-38},{-8,-38}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction1.products[1], boundary_fore1.rear) annotation (Line(
          points={{12,-38},{22,-38}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Test for the undirected flow resistance.</u>
<u><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
    end TestFlowResistance;

    model SimpleReaction "The simple chemical reaction A<->B with equilibrium B/A = 2"
      import Chemical;
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(
        useRear=false,
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

      Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(K))},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Chemical.Boundaries.Substance B(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{42,-8},{62,12}})));

      inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{58,64},{78,84}})));
    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-48,-8},{-48,-92},{60,-92},{60,-98}},
          color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{46,-8},{46,-92},{60,-92},{60,-98}},
                                         color={127,127,0}));
      connect(A.fore, reaction2_1.substrates[1]) annotation (Line(
          points={{-32,2},{-10,2}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], B.rear) annotation (Line(
          points={{10,2},{42,2}},
          color={158,66,200},
          thickness=0.5));
      annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end SimpleReaction;

    model SimpleReaction2 "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
      import Chemical;
      import Chemical;
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(
        useFore=true,
        useSolution=true,
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
      Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=2, DfG=-R*T_25degC*log(Kx))},
        nS=2,
        nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Boundaries.Substance B(
        useFore=true,
        useSolution=true,
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Boundaries.Substance C(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,4},{-30,-90},{60,-90},{60,-98}},
          color={127,127,0}));
      connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                 color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
            {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

      connect(A.fore, reaction2_1.substrates[1])
        annotation (Line(
          points={{-14,14},{-4,14},{-4,1.75},{4,1.75}},
          color={158,66,200},
          thickness=0.5));
      connect(B.fore, reaction2_1.substrates[2])
        annotation (Line(
          points={{-14,-14},{-8,-14},{-8,2.25},{4,2.25}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], C.rear) annotation (Line(
          points={{24,2},{48,2}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end SimpleReaction2;

    model SimpleReaction22 "The simple chemical reaction A+B<->C+D with equilibrium [C]*[D]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Chemical.Processes.Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx)),
            Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx))},
        nS=2,
        nP=2) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Boundaries.Substance B(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Boundaries.Substance C(amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

      Chemical.Boundaries.Substance D(amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{44,-34},{64,-14}})));
      inner DropOfCommons dropOfCommons(L=1e-3) annotation (Placement(transformation(extent={{52,56},{72,76}})));
    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,2},{-30,-90},{60,-90},{60,-98}},
          color={127,127,0}));
      connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                 color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
            {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

      connect(B.outlet, reaction2_1.substrates[1]) annotation (Line(
          points={{-14,-14},{-4,-14},{-4,1.75},{4,1.75}},
          color={158,66,200},
          thickness=0.5));
      connect(A.outlet, reaction2_1.substrates[2]) annotation (Line(
          points={{-14,12},{-4,12},{-4,2.25},{4,2.25}},
          color={158,66,200},
          thickness=0.5));
      connect(D.solution, solution.solution) annotation (Line(points={{48,-34},{60,-34},{60,-98}}, color={127,127,0}));
      connect(reaction2_1.products[1], C.inlet) annotation (Line(
          points={{24,1.75},{36,1.75},{36,2},{48,2}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[2], D.inlet) annotation (Line(
          points={{24,2.25},{34,2.25},{34,-24},{44,-24}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
    end SimpleReaction22;

    model ExothermicReaction "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;

       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

      Chemical.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
        annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
      Chemical.Boundaries.Substance A(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Reaction reaction2_2(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Chemical.Boundaries.Substance B(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

      Chemical.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,0},{98,94}})));
      Chemical.Boundaries.Substance A1(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Chemical.Boundaries.Substance B1(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{20,40},{40,60}})));

      //  Modelica.SIunits.HeatFlowRate q
      //    "Heat flow to environment to reach constant temperature";
      Modelica.Units.SI.Temperature t
        "Temperature if the solution is ideally thermal isolated from environment";
      Chemical.Boundaries.Substance H2O(substanceData=Chemical.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{20,4},{40,24}})));
      Chemical.Boundaries.Substance H2O1(substanceData=Chemical.Substances.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{20,-94},{40,-74}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,4},{-8,24}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed2
        annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
      inner Modelica.Fluid.System system(T_ambient=298.15)
        annotation (Placement(transformation(extent={{56,64},{76,84}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
            298.15)
        annotation (Placement(transformation(extent={{-88,26},{-68,46}})));
    equation
      //  q = fixedTemperature.port.Q_flow;
      t = thermal_isolated_solution.solution.T;

      connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{24,-60},{24,-64},{58.4,-64},{58.4,-99.06}},
          color={127,127,0}));
      connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
            points={{-36,-60},{-36,-64},{58.4,-64},{58.4,-99.06}},
                                                             color={127,127,0}));
      connect(B1.solution, solution_at_constant_temperature.solution) annotation (
          Line(
          points={{24,40},{24,34},{58.4,34},{58.4,0.94}},
          color={127,127,0}));
      connect(A1.solution, solution_at_constant_temperature.solution) annotation (
          Line(points={{-36,40},{-36,34},{58.4,34},{58.4,0.94}},
                                                          color={127,127,0}));
    connect(solution_at_constant_temperature.solution, H2O.solution)
      annotation (Line(
        points={{58.4,0.94},{24,0.94},{24,4}},
        color={127,127,0}));
    connect(thermal_isolated_solution.solution, H2O1.solution) annotation (Line(
        points={{58.4,-99.06},{24,-99.06},{24,-94}},
        color={127,127,0}));
    connect(solution_at_constant_temperature.bottom, fixed1.flange) annotation (
       Line(
        points={{-1,-0.94},{0,-0.94},{0,14},{-18,14}},
        color={0,127,0}));
    connect(thermal_isolated_solution.bottom, fixed2.flange) annotation (Line(
        points={{-1,-100.94},{-1,-86},{-16,-86}},
        color={0,127,0}));
      connect(solution_at_constant_temperature.heatPort, fixedTemperature.port)
        annotation (Line(points={{-60.4,-0.94},{-60.4,36},{-68,36}}, color={191,0,
              0}));
      connect(A1.fore, reaction2_1.substrates[1]) annotation (Line(
          points={{-20,50},{-8,50}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], B1.rear) annotation (Line(
          points={{12,50},{20,50}},
          color={158,66,200},
          thickness=0.5));
      connect(A.fore, reaction2_2.substrates[1]) annotation (Line(
          points={{-20,-50},{-8,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_2.products[1], B.rear) annotation (Line(
          points={{12,-50},{20,-50}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}})));
    end ExothermicReaction;

    model ConductionElement "Test for ConductionElement"
      extends Modelica.Icons.Example;

      replaceable package Medium =
          Chemical.Media.myMedia.Incompressible.Examples.Glycol47                          constrainedby
        Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medium Model"
        annotation(choicesAllMatching=true, Documentation(info="<html>
<u>
Medium Model for the test. Be aware that the Component is mainly
meant for liquids with low compressablility.
</u>
</html>"));

      Chemical.Processes.ConductionElement conductionElement(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.T,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,60},{10,80}})));

      Chemical.Boundaries.BoundaryFore boundary_fore(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
      Modelica.Blocks.Sources.Step step(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=373.15)
        annotation (Placement(transformation(extent={{80,50},{60,70}})));
      inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{60,-70},{80,-50}})));
      Chemical.Processes.ConductionElement conductionElement1(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.h,
        T_0=263.15,
        h_0=1000,
        neglectChemicalPotentialChanges=false) annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Chemical.Boundaries.BoundaryFore boundary_fore1(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,30},{40,50}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=2e2,
        duration=0.001,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,30},{-50,50}})));
      Chemical.Processes.ConductionElement conductionElement2(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        enforce_global_energy_conservation=true,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.fore,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,0},{10,20}})));

      Chemical.Boundaries.BoundaryFore boundary_fore2(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,0},{40,20}})));
      Chemical.Boundaries.BoundaryRear boundary_rear2(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
      Modelica.Blocks.Sources.Step step1(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Chemical.Processes.ConductionElement conductionElement3(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.rear,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

      Chemical.Boundaries.BoundaryFore boundary_fore3(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-30},{40,-10}})));
      Chemical.Boundaries.BoundaryRear boundary_rear3(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      Modelica.Blocks.Sources.Step step2(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
      Chemical.Processes.ConductionElement conductionElement4(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=0,
        V(displayUnit="l") = 0.001,
        enforce_global_energy_conservation=true,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.port,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Chemical.Boundaries.BoundaryFore boundary_fore4(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
      Chemical.Boundaries.BoundaryRear boundary_rear4(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Modelica.Blocks.Sources.Step step3(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
      Chemical.Processes.ConductionElement conductionElement5(
        redeclare package stateOfMatter = stateOfMatter,
        L=0.2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        n_flow_0=-1,
        V(displayUnit="l") = 0.001,
        A=35,
        U=500,
        init=Chemical.Processes.Internal.InitializationMethodsCondElement.port,
        T_0=263.15) annotation (Placement(transformation(extent={{-10,-90},{10,-70}})));

      Chemical.Boundaries.BoundaryFore boundary_fore5(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=328.15,
        u0_par=100000) annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
      Chemical.Boundaries.BoundaryRear boundary_rear5(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=288.15) annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
      Modelica.Blocks.Sources.Step step4(
        height=2e2,
        offset=0.999e5,
        startTime=0.33)
        annotation (Placement(transformation(extent={{-70,-90},{-50,-70}})));
    equation
      connect(step.y, boundary_rear.u0_var) annotation (Line(points={{-49,70},{-40,70},{-40,76},{-32,76}},
                                              color={0,0,127}));
      connect(conductionElement.rear, boundary_rear.fore) annotation (Line(
          points={{-10,70},{-20,70}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement.fore, boundary_fore.rear) annotation (Line(
          points={{10,70},{20,70}},
          color={158,66,200},
          thickness=0.5));
      connect(ramp.y, boundary_rear1.u0_var) annotation (Line(points={{-49,40},{-40,40},{-40,46},{-32,46}},
                                                                                            color={0,0,127}));
      connect(conductionElement1.rear, boundary_rear1.fore)
        annotation (Line(
          points={{-10,40},{-20,40}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement1.fore, boundary_fore1.rear) annotation (Line(
          points={{10,40},{20,40}},
          color={158,66,200},
          thickness=0.5));
      connect(step1.y, boundary_rear2.u0_var) annotation (Line(points={{-49,10},{-40,10},{-40,16},{-32,16}},
                                                                                           color={0,0,127}));
      connect(conductionElement2.rear, boundary_rear2.fore) annotation (Line(
          points={{-10,10},{-20,10}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement2.fore, boundary_fore2.rear) annotation (Line(
          points={{10,10},{20,10}},
          color={158,66,200},
          thickness=0.5));
      connect(step2.y, boundary_rear3.u0_var) annotation (Line(points={{-49,-20},{-40,-20},{-40,-14},{-32,-14}},
                                                                                             color={0,0,127}));
      connect(conductionElement3.rear, boundary_rear3.fore)
        annotation (Line(
          points={{-10,-20},{-20,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement3.fore, boundary_fore3.rear) annotation (Line(
          points={{10,-20},{20,-20}},
          color={158,66,200},
          thickness=0.5));
      connect(step3.y, boundary_rear4.u0_var) annotation (Line(points={{-49,-50},{-40,-50},{-40,-44},{-32,-44}},
                                                                                             color={0,0,127}));
      connect(conductionElement4.rear, boundary_rear4.fore)
        annotation (Line(
          points={{-10,-50},{-20,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement4.fore, boundary_fore4.rear) annotation (Line(
          points={{10,-50},{20,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(step4.y, boundary_rear5.u0_var) annotation (Line(points={{-49,-80},{-40,-80},{-40,-74},{-32,-74}},
                                                                                             color={0,0,127}));
      connect(conductionElement5.rear, boundary_rear5.fore)
        annotation (Line(
          points={{-10,-80},{-20,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(conductionElement5.fore, boundary_fore5.rear) annotation (Line(
          points={{10,-80},{20,-80}},
          color={158,66,200},
          thickness=0.5));
      connect(fixedTemperature.port, conductionElement1.heatPort) annotation (Line(points={{60,60},{40,60},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      connect(conductionElement.heatPort, conductionElement1.heatPort)
        annotation (Line(points={{0,79.8},{0,88},{40,88},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      connect(conductionElement2.heatPort, conductionElement1.heatPort)
        annotation (Line(points={{0,19.8},{0,26},{40,26},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      connect(conductionElement3.heatPort, conductionElement1.heatPort)
        annotation (Line(points={{0,-10.2},{0,-4},{40,-4},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      connect(conductionElement4.heatPort, conductionElement1.heatPort)
        annotation (Line(points={{0,-40.2},{0,-34},{40,-34},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      connect(conductionElement5.heatPort, conductionElement1.heatPort)
        annotation (Line(points={{0,-70.2},{0,-66},{40,-66},{40,56},{0,56},{0,49.8}}, color={191,0,0}));
      annotation (experiment(StopTime=1, Tolerance=1e-6, Interval=0.001),
      Documentation(info="<html>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
    end ConductionElement;

    model TransportDelay "Test for transport delay"
      extends Modelica.Icons.Example;

      replaceable package Medium = Chemical.Media.myMedia.Air.DryAirNasa constrainedby Chemical.Media.myMedia.Interfaces.PartialMedium
                                                             "Medium Model"
        annotation (Documentation(info="<html>
<u>
Medium model for the test. Can be anything.
</u>
</html>"));

      Chemical.Processes.TransportDelay transportDelay(
        redeclare package stateOfMatter = stateOfMatter,
        l=100,
        r(displayUnit="mm") = 0.015)
        annotation (Placement(transformation(extent={{30,30},{50,50}})));
      Chemical.Processes.FlowResistance flowResistance(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=100,
        l(displayUnit="mm") = 0.008,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
            (
          k=1e4))
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Chemical.Boundaries.Sink sink(redeclare package stateOfMatter =
            stateOfMatter,
          u0_par=100000)
        annotation (Placement(transformation(extent={{70,30},{90,50}})));
      Chemical.Boundaries.Source source(
        redeclare package stateOfMatter = stateOfMatter,
        temperatureFromInput=true,
        potentialFromInput=true)
        annotation (Placement(transformation(extent={{-48,30},{-28,50}})));
      Modelica.Blocks.Sources.Trapezoid
                                   trapezoid(
        amplitude=-6e3,
        rising=0.2,
        width=0.2,
        falling=0.2,
        period=0.8,
        offset=1.04e5,
        startTime=0.2)
        annotation (Placement(transformation(extent={{-94,10},{-74,30}})));
      Modelica.Blocks.Sources.Ramp ramp1(
        height=20,
        duration=1.5,
        offset=273,
        startTime=0.3)
        annotation (Placement(transformation(extent={{-94,50},{-74,70}})));
      inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{60,-80},{80,-60}})));

      Chemical.Processes.TransportDelay transportDelay1(
        redeclare package stateOfMatter = stateOfMatter,
        l=100,
        r(displayUnit="mm") = 0.015) annotation (Placement(transformation(extent={{30,-50},{50,-30}})));
      FlowResistance flowResistance1(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=100,
        l(displayUnit="mm") = 0.008,
        redeclare function pLoss =
            Chemical.Processes.Internal.FlowResistance.linearQuadraticChemicalPotentialLoss
            (
          k=1e4))
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
      Chemical.Boundaries.BoundaryFore boundary_fore(redeclare package stateOfMatter = stateOfMatter, u0_par=100000)
        annotation (Placement(transformation(extent={{70,-50},{90,-30}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
        redeclare package stateOfMatter = stateOfMatter,
        temperatureFromInput=true,
        potentialFromInput=true) annotation (Placement(transformation(extent={{-48,-50},{-28,-30}})));
      Modelica.Blocks.Sources.Trapezoid trapezoid1(
        amplitude=-6e3,
        rising=0.2,
        width=0.2,
        falling=0.2,
        period=0.8,
        offset=1.04e5,
        startTime=0.2)
        annotation (Placement(transformation(extent={{-94,-70},{-74,-50}})));
      Modelica.Blocks.Sources.Ramp ramp2(
        height=20,
        duration=1.5,
        offset=273,
        startTime=0.3)
        annotation (Placement(transformation(extent={{-94,-30},{-74,-10}})));
    equation
      connect(ramp1.y, source.T0_var) annotation (Line(points={{-73,60},{-60,60},{-60,40},{-40,40}},
                               color={0,0,127}));
      connect(trapezoid.y, source.u0_var) annotation (Line(points={{-73,20},{-60,20},{-60,46},{-40,46}},
                                      color={0,0,127}));
      connect(flowResistance.inlet, source.outlet) annotation (Line(
          points={{-10,40},{-28,40}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance.outlet, transportDelay.inlet) annotation (Line(
          points={{10,40},{30,40}},
          color={158,66,200},
          thickness=0.5));
      connect(transportDelay.outlet, sink.inlet) annotation (Line(
          points={{50,40},{70,40}},
          color={158,66,200},
          thickness=0.5));
      connect(ramp2.y, boundary_rear.T0_var) annotation (Line(points={{-73,-20},{-60,-20},{-60,-40},{-40,-40}},
                                              color={0,0,127}));
      connect(trapezoid1.y, boundary_rear.u0_var) annotation (Line(points={{-73,-60},{-60,-60},{-60,-34},{-40,-34}},
                                                   color={0,0,127}));
      connect(flowResistance1.rear, boundary_rear.fore) annotation (Line(
          points={{-10,-40},{-28,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance1.fore, transportDelay1.rear) annotation (Line(
          points={{10,-40},{30,-40}},
          color={158,66,200},
          thickness=0.5));
      connect(transportDelay1.fore, boundary_fore.rear) annotation (Line(
          points={{50,-40},{62,-40},{62,-40},{70,-40}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        experiment(
          StopTime=2.5,
          Tolerance=1e-6,
          Interval=0.0025,
          __Dymola_Algorithm="Dassl"),
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Test for use Transport Delay to delay the thermodynamic state of the flow. </u>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
    end TransportDelay;

    model EnzymeKinetics "Basic enzyme kinetics"
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;
      import Chemical;
      extends Modelica.Icons.Example;

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Chemical.Boundaries.Substance P(
        useRear=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{72,-12},{92,8}})));

      Chemical.Boundaries.Substance S(
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=100) annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

      parameter Modelica.Units.SI.AmountOfSubstance tE=1
        "Total amount of enzyme";
         parameter Real k_cat(
        unit="mol/s",
        displayUnit="mol/min")=1
        "Forward rate of second reaction";
      constant Modelica.Units.SI.Concentration Km=0.1
        "Michaelis constant = substrate concentration at rate of half Vmax";

      parameter Modelica.Units.SI.MolarFlowRate Vmax=1e-5*k_cat
        "Maximal molar flow";

      Chemical.Boundaries.Substance ES(
        useRear=true,
        useFore=true,
        useSolution=true,
        initAmount=Chemical.Utilities.Types.InitializationMethods.state,
        use_mass_start=false,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
      Chemical.Boundaries.Substance E(
        useRear=true,
        useFore=true,
        useSolution=true,
        use_mass_start=false,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{12,36},{-8,56}})));
      Reaction                    chemicalReaction(
        k_forward=1,
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*log(2/Km))},
              nS=2,
        nP=1) annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));

      Processes.ForwardReaction chemicalReaction1(
        k_forward=k_cat,
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*50),
            Chemical.Interfaces.Incompressible.SubstanceDataParameters()},
        nS=1,
        nP=2) annotation (Placement(transformation(extent={{24,-8},{44,12}})));

      Chemical.Boundaries.Substance liquidWater(
        useSolution=true,
        substanceData=Chemical.Substances.Water_liquid(),
        use_mass_start=true,
        mass_start=1) annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
      inner DropOfCommons dropOfCommons      annotation (Placement(transformation(extent={{68,70},{88,90}})));
    equation
      //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);
      connect(E.solution, solution.solution) annotation (Line(
          points={{8,36},{-8,36},{-8,-98},{60,-98}},
          color={127,127,0}));
      connect(ES.solution, solution.solution)
        annotation (Line(points={{-4,-10},{-4,-98},{60,-98}},         color={127,127,0}));

      connect(S.solution, solution.solution) annotation (Line(
          points={{-88,-14},{-88,-56},{-8,-56},{-8,-98},{60,-98}},
          color={127,127,0}));
      connect(P.solution, solution.solution) annotation (Line(
          points={{76,-12},{76,-98},{60,-98}},
          color={127,127,0}));
      connect(liquidWater.solution, solution.solution) annotation (Line(points={{
              46,-80},{46,-98},{60,-98}}, color={127,127,0}));
      connect(chemicalReaction.products[1], ES.rear) annotation (Line(
          points={{-22,0},{-8,0}},
          color={158,66,200},
          thickness=0.5));
      connect(ES.fore, chemicalReaction1.substrates[1])
        annotation (Line(
          points={{12,0},{18,0},{18,2},{24,2}},
          color={158,66,200},
          thickness=0.5));
      connect(S.fore, chemicalReaction.substrates[1])
        annotation (Line(
          points={{-72,-4},{-52,-4},{-52,-2},{-42,-2},{-42,-0.25}},
          color={158,66,200},
          thickness=0.5));
      connect(E.fore, chemicalReaction.substrates[2])
        annotation (Line(
          points={{-8,46},{-52,46},{-52,0.25},{-42,0.25}},
          color={158,66,200},
          thickness=0.5));
      connect(chemicalReaction1.products[1], P.rear)
        annotation (Line(
          points={{44,1.75},{58,1.75},{58,-2},{72,-2}},
          color={158,66,200},
          thickness=0.5));
      connect(E.rear, chemicalReaction1.products[2])
        annotation (Line(
          points={{12,46},{48,46},{48,48},{56,48},{56,2.25},{44,2.25}},
          color={158,66,200},
          thickness=0.5));
          annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Be carefull, the assumption for Michaelis-Menton are very strong: </p>
<p>The substrate must be in sufficiently high concentration and the product must be in very low concentration to reach almost all enzyme in enzyme-substrate complex all time. ([S] &gt;&gt; Km) &amp;&amp; ([P] &lt;&lt; K2)</p>
<p><br>To recalculate the enzyme kinetics from Michaelis-Menton parameters Km, tE a k_cat is selected the same half-rate of the reaction defined as:</p>
<p>E = ES = tE/2 .. the amount of free enzyme is the same as the amount of enzyme-substrate complexes</p>
<p>S = Km .. the amount of substrate is Km</p>
<p>r = Vmax/2 = tE*k_cat / 2 .. the rate of reaction is the half of maximal rate</p>
<p><br>Conversions of molar concentration to mole fraction (MM is molar mass of the solvent in solution -&gt; 55.508 kg/mol for water):</p>
<p>x(Km) = Km/MM</p>
<p>x(tE) = tE/MM</p>
<p>xS = S/MM = Km/MM</p>
<p><br>The new kinetics of the system defined as:</p>
<p>uS&deg; = DfG(S) = 0</p>
<p>uE&deg; = DfG(E) = 0</p>
<p>uES&deg; = <b>DfG(ES) = DfG(S) + DfG(E) - R*T*ln(2/x(Km))</b></p>
<p>from dissociation coeficient of the frist reaction 2/x(Km) = xSE/(xS*xE) = exp((uE&deg; + uS&deg; - uES&deg;)/(RT))</p>
<p>uP&deg; = DfG(P) </p>
<p><br>r = Vmax/2</p>
<p>r = -kC1 * (uES&deg; - uE&deg; - uS&deg; + R*T*ln(xES/(xE*xS) ) = -kC1 * (-R*T*ln(2/x(Km)) + R*T*ln(xS) ) = kC1 * R * T * ln(2)</p>
<p>because xES=xE this time</p>
<p>r = -kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) = -kC2 * (DfG(P) - uES&deg; + R*T*ln(xP) ) = kC2 * (-DfG(P) - R * T * ln(2))</p>
<h4>kC1 = (Vmax/2) / (R * T * ln(2))</h4>
<h4>kC2 = (Vmax/2) / ( -DfG(P) - R * T * ln(2) ) </h4>
<p><br>For example in case of C=AmountOfSolution/(Tau*ActivationPotential) we can rewrite C to ActivationPotential (Be carefull: this energy is not the same as in <a href=\"http://en.wikipedia.org/wiki/Arrhenius_equation\">Arrhenius equation</a> or in Transition State Theory):</p>
<p>ActivationPotential1 = AmountOfSolution/(Tau*(Vmax/2)) * R * T * ln(2) </p>
<p>ActivationPotential2 = AmountOfSolution/(Tau*(Vmax/2)) * ( -DfG(P) - R * T * ln(2) ) </p>
<p><br>where</p>
<p>AmountOfSolution = MM = 55.508 (for water)</p>
<p>Tau = 1 s (just to be physical unit correct)</p>
<p>DfG(P) = -R*T*50 is Gibbs energy of formation of product (setting negative enough makes second reaction almost irreversible)</p>
<h4>The maximum of the new enzyme kinetics</h4>
<p>The enzymatic rate must have a maximum near of Vmax. </p>
<p>The new maximum is a litle higher: Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) ), for example if -uP&deg;/RT = 50, the new maximum is around 1.014*Vmax, where Vmax is the maximum of Michaelis Menten.</p>
<p>The proof:</p>
<p>We want to sutisfied the following inequality:</p>
<p>-kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) ?=&lt;? Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) )</p>
<p><br>(Vmax/2) * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) / ( - uP&deg; - R * T * ln(2) ) ?=&lt;? Vmax*(1 + R*T*ln(2) / ( -uP&deg; - R*T*ln(2)) )</p>
<p>(uP&deg; +<b> </b>R*T*ln(2/x(Km)) + R*T*ln(xP*xE/xES) ) ?=&lt;? 2*( - uP&deg; - R * T * ln(2) ) + 2*R*T*ln(2)</p>
<p>R*T*ln(xP*xE/xES) ?=&lt;? - uP&deg; - R*T*ln(2/x(Km)) </p>
<p>xP*xE/xES ?=&lt;? exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>The equality is the equation of the equilibrium: xP*xE/xES = exp((- uP&deg; - uE&deg; + uES&deg; )/(R*T)) = exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>If the equilibrium of the reaction is reached only by forward rate then xP*xE/xES must be less than the dissociation constant.</p>
</html>"),
        experiment(StopTime=100000, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end EnzymeKinetics;
    annotation (Documentation(info="<html>
<u>Tests for top level components of the undirected chemical simulation package.</u>
</html>"));
  end Tests;

  package Internal "Internals package for Processes"
    extends Modelica.Icons.InternalPackage;

    partial model PartialConductionElement "Partial volume with quasisationary mass and heatport and undetermined heat transfer coefficient"
      extends Chemical.Interfaces.SISO                 (
                                    final cliu_u_out=false);

      parameter Modelica.Units.SI.Volume V=1 "Volume of the element";
      parameter Chemical.Processes.Internal.InitializationMethodsCondElement init=Chemical.Processes.Internal.InitializationMethodsCondElement.rear
        "Initialization method for h" annotation (Dialog(tab="Initialization", group="Enthalpy"));
      parameter Modelica.Units.SI.Temperature T_0=Medium.T_default "Initial Temperature" annotation (Dialog(
          tab="Initialization",
          group="Enthalpy",
          enable=(init == Chemical.Processes.Internal.InitializationMethodsCondElement.T)));
      parameter Modelica.Units.SI.MolarEnthalpy h_0=Medium.h_default "Initial molar enthalpy" annotation (Dialog(
          tab="Initialization",
          group="Enthalpy",
          enable=(init == Chemical.Processes.Internal.InitializationMethodsCondElement.h)));
      parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimal density" annotation (Dialog(tab="Advanced"));
      parameter Boolean neglectChemicalPotentialChanges = true "Neglect potential changes in energy equation"
        annotation(Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg
        "Regularization massflow to switch between positive- and negative-massflow model" annotation (Dialog(tab="Advanced"));
      parameter Boolean enforce_global_energy_conservation = false "If true, exact global energy conservation is enforced by feeding back all energy stored locally back in the system"
        annotation(Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.Time T_e=100 "Factor for feeding back energy."
        annotation (Dialog(tab="Advanced", enable=enforce_global_energy_conservation));

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(Q_flow=Q_flow, T=T_heatPort)
        annotation (Placement(transformation(extent={{-10,88},{10,108}})));

      Modelica.Units.SI.MolarEnthalpy h(start=Medium.h_default, stateSelect=StateSelect.prefer);

      Medium.ThermodynamicState state = Chemical.Interfaces.SubstanceState(u=u,h= h);
      Modelica.Units.SI.Temperature T=Medium.temperature(state);
      Modelica.Units.SI.ThermalConductance k "Thermal conductance heatport->fluid";

      Modelica.Units.SI.Energy deltaE_system(start=0, fixed=true) "Energy difference between n_flow*(h_in-h_out) and Q_flow";

      Modelica.Units.SI.Mass M;

    protected
      Modelica.Units.SI.ChemicalPotential u=Chemical.Utilities.Internal.regStep(
                n_flow,
                u_rear_in,
                u_fore_in,
                n_flow_reg);

      //h_in only in rhs of ODE--> h still smooth, better results at low massflow than using regStep
      Modelica.Units.SI.MolarEnthalpy h_in=if n_flow >= 0 then h_rear_in else h_fore_in;

      Modelica.Units.SI.Density rho=max(rho_min, Medium.density(state));
      Modelica.Units.SI.MolarEnthalpy h_out=state.h;

      Modelica.Units.SI.Temperature T_heatPort;
      Modelica.Units.SI.HeatFlowRate Q_flow;

    initial equation
      if init == Internal.InitializationMethodsCondElement.T then
         h = Medium.setState_pTX(u, T_0, Xi);
      elseif init == Internal.InitializationMethodsCondElement.h then
         h = h_0;
      elseif init == InitializationMethodsCondElement.port then
         h = h_in;
      elseif init == Internal.InitializationMethodsCondElement.fore then
         h = h_fore_in;
       else
         h = h_rear_in;
       end if;

    equation
      //mass is assumed to be quasisationary-> this violates the conservation of mass since n_flow_in = -n_flow_out. see documentation/information
      M = V*rho;

      // Mass is assumed to be (quasi-)stationary-> assumption: der(M)=0.
      // V=const
      // lhs: der(U) = der(H - pV) = der(M)*h +  M*der(h) - V*der(u) = M*der(h)
      // one can not simply add der(M)*h to the left side, since the right side is also assuming sationary mass (n_flow_in=-n_flow_out)
      // on RHS: Q_flow + n_flow_in*h_in + n_flow_out*h = Q_flow + n_flow_in*h_in + (n_flow_in-n_flow_in+ n_flow_out)*h = Q_flow + n_flow_in*(h_in-h) + der(M)*h = Q_flow + n_flow*(h_in-h)
      // or:     Q_flow + n_flow_in*h_in + n_flow_out*h = Q_flow + (n_flow_in+n_flow_out-n_flow_out)*h_in + n_flow_out*h = Q_flow + n_flow_out*(h-h_in) + der(M)*h_in = Q_flow + n_flow*(h_in - h)
      // so there is a resulting term on the RHS with der(M)*(h(k) + h_out(1-k)) with k in [0,1]; that also would have to be accounted for when der(M) not 0 on LHS
      // to improve the equation u_in and M would have to be filtered and a value for k would have to be found, all that resulted in only moderate improvements
      // therefore energy that is accumulated with system border around the component (deltaE_system) is slowly fed back into the system when enforce_global_energy_conservation is true
      //(then energy stored in the system e.g. by evaproation/condension, while still being visible short-term, neglegted in logterm)
      if not neglectChemicalPotentialChanges then
        M*der(h) = Q_flow + abs(n_flow)*(h_in - h_out) + V*der(u) + (if enforce_global_energy_conservation then deltaE_system/T_e else 0);
      else
        //neglegt V*der(u), since u might not be smooth -> to notice a difference der(u) must be about 1e7 Pa/s. see documentation/information
        M*der(h) = Q_flow + abs(n_flow)*(h_in - h_out) + (if enforce_global_energy_conservation then deltaE_system/T_e else 0);
      end if;
      if enforce_global_energy_conservation then
        der(deltaE_system) = Q_flow + abs(n_flow)*(h_in - h_out);
      else
        deltaE_system = 0;
      end if;

      Q_flow = k*(T_heatPort - T);

      //forwards model
      du_fore = 0;
      h_fore_out = h;
      Xi_fore_out = Xi_rear_in;

      //rearwards model
      du_rear = 0;
      h_rear_out = h;
      Xi_rear_out = Xi_fore_in;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
           Line(
             points={{-100,0},{100,0}},
             thickness=0.5,
             color={158,66,200}),
           Ellipse(
             extent={{-70,-70},{70,70}},
             lineColor={158,66,200},
             lineThickness=0.5,
             fillColor={194,138,221},
             fillPattern=FillPattern.Solid,
             pattern=LinePattern.Solid),
           Line(
             points={{-50,-30},{50,-30}},
             color={238,46,47}),
           Line(
             points={{-50,-15},{50,-15}},
             color={238,46,47}),
           Line(
             points={{-50,0},{50,0}},
             color={238,46,47}),
           Line(
             points={{-50,15},{50,15}},
             color={238,46,47}),
           Line(
             points={{-50,30},{50,30}},
             color={238,46,47}),
           Line(
             points={{0,100},{0,-30}},
             color={238,46,47})}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>
Undirected implementation of the Conduction Element.
</u>
<u>
This model is an element with a fixed volume (fig. 1). The mass in the volume is 
assumed quasi-stationary (statically computed with volume and density), and the 
fore massflow is coupled to the rear massflow. <strong>Because of this the 
ConductionElement cannot be used as a loop breaker</strong>.
The advantage is that multiple ConductionElements can be put behind each other 
without worrying about oscillations or fast eigenvalues between their masses. 
The ConductionElement implements equations for conservation of mass and energy 
for the fluid mass contained within it.
</u>
<u>
Different to the unidirectional ConductionElement, the model for forward massflow
(see fig. 1) is valid for both flow directions.
</u>
<u>
Additionally the undirected Conduction Element offers more initialization methods 
for h. Additionally to initialize T or h by a parameter, one can choose to 
initialize with the incoming enthalpy from either one of the two ports, or use the 
correct one, depending on the massflow (option &apos;port&apos;). The last option 
can lead to large nonlinear initialization problems, we advice to choose the port 
to initialize h from if known in advance (options &apos;rear&apos; or &apos;fore&apos;).
</u>
<u>
The ConductionElement makes different assumptions:
</u>
<ul>
  <li>
    Quasistationary Mass:
    <br>
    n_flow_rear = - n_flow_fore &amp; M = rho * V
    (this assumption violates the conservation of mass for changing densities,
    since the mass in the element can change although inflow and outflow are the
    same)
    <br>
    der(H) = der(M*h) = M*der(h)
    (This assumption violates the conservation of energy for changing densities,
    since then the mass M of fluid in the element is no longer constant)
  </li>
  <li>
    Neglection of der(u) in the energy equation
    <br>
    V*der(u) = 0
    (this assumption violates the conservation of energy for changing potentials.
    For a noticeable difference in the testcase the der(u) must be in the
    order of 1e5 Pa/s).
    <br>
    This assumption can be turned off by setting neglectChemicalPotentialChanges=false
    (true by default) in the Advanced tab. <strong>This option requires the 
    fore and rear input potentials to be smooth.</strong>
  </li>
</ul>
<u>
Due to these assumptions minor violations in the global energy conservation can 
occur.
With the flag enforce_global_energy_conservation in the &quot;Advanced&quot; tab 
is set true (Default: false), long-term energy storage in the ConductionElement 
is sacrificed to hold global energy conservation.
</u>
<div>
<img src=\"modelica://Chemical/Resources/Doku/Chemical.Processes.ConductionElement_positive.png\"/>
</div>
<u>
fig. 1: positive massflow model
</u>
</html>"));
    end PartialConductionElement;

    type InitializationMethodsCondElement = enumeration(
        T "Temperature T_0",
        h "molar enthalpy h_0",
        fore "input state from fore",
        rear "input state from rear",
        port "regularized input state from fore or rear, depending on massflow") "Choices for initialization of state h of undirected ConductionElement"
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>
Choices for initialization of a state h.
</u>
</html>"));
    partial model PartialReactionWithSubstanceData "Chemical Reaction"
      extends Chemical.Processes.Internal.PartialReaction;

      parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP]
       annotation (choicesAllMatching = true);

    equation

      products.definition = productsSubstanceData;
      if nS>0 then
        products.solution = fill(substrates[1].solution,nP);
      end if;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
<p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
<p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"modelica://ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
<td><p>dissociation constant</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
<td><p>molar Gibb&apos;s energy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>H = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>H<sub>i</sub>) </p></td>
<td><p>molar enthalpy of the reaction</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>S = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>S<sub>i</sub>) = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>r</sub>&omega;) </p></td>
<td><p>molar entropy of the reaction</p></td>
</tr>
</table>
<h4><span style=\"color:#008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stochiometric coefficients of i-th substance</p></td>
</tr>
<tr>
<td><p>K</p></td>
<td><p>dissociation constant (activity based)</p></td>
</tr>
<tr>
<td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
<td><p>activity of the substance A</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub></p></td>
<td><p>activity coefficient of the substance A</p></td>
</tr>
<tr>
<td><p>x<sub>i</sub></p></td>
<td><p>mole fraction of the substance A</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
<td><p>molar enthalpy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
<td><p>molar Gibbs energy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
<td><p>molar entropy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by reaction</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
    end PartialReactionWithSubstanceData;

    partial model PartialReaction "Chemical Reaction"
      import Chemical;
      import Chemical.Utilities.Types.InitializationMethods;

      replaceable package stateOfMatter = Chemical.Interfaces.Incompressible constrainedby Chemical.Interfaces.StateOfMatter
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

      parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
        annotation(Dialog(tab="Advanced"));
      parameter InitializationMethods initN_flow =Chemical.Utilities.Types.InitializationMethods.none  "Initialization method for n_flow"
        annotation(Dialog(tab= "Initialization", group="Molar flow"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
        annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
      parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
        annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

      parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
      parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
        annotation(Dialog(tab="Advanced"));

      parameter Integer nS=0 "Number of substrate types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

      parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
        "Stoichiometric reaction coefficient for substrates"
        annotation (HideResult=true);

      parameter Integer nP=0 "Number of product types"
        annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

      parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
        "Stoichiometric reaction coefficients for products"
        annotation (HideResult=true);

      Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Reaction molar flow rate";

      Chemical.Interfaces.Rear substrates[nS](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-100,0})));

      Chemical.Interfaces.Fore products[nP](redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0}), iconTransformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={100,0})));

      Modelica.Units.SI.MolarEnthalpy h_fore_mix, h_rear_mix;

      Real duRT_fore, duRT_rear, du_fore, du_rear, dr, Sx_fore,Px_rear,Kx;

      Modelica.Units.SI.ChemicalPotential uPure_substrates[nS];
      Modelica.Units.SI.ChemicalPotential uPure_products[nP];
    protected
      outer DropOfCommons dropOfCommons;

    public
        //debug only:
      Real Px_fore,Sx_rear;

    initial equation
      if initN_flow == InitializationMethods.state then
        rr = n_flow_0;
      elseif initN_flow == InitializationMethods.derivative then
        der(rr) = n_acceleration_0;
      elseif initN_flow == InitializationMethods.steadyState then
        der(rr) = 0;
      end if;

    equation

      du_fore = (s * substrates.state_forwards.u) - (p * products.state_forwards.u);
      du_rear = (s * substrates.state_rearwards.u) - (p * products.state_rearwards.u);

      duRT_fore = ((s * (substrates.state_forwards.u ./ (Modelica.Constants.R*substrates.solution.T))) - (p * (products.state_forwards.u ./ (Modelica.Constants.R*products.solution.T))));
      duRT_rear = ((s * (substrates.state_rearwards.u ./ (Modelica.Constants.R*substrates.solution.T))) - (p * (products.state_rearwards.u ./ (Modelica.Constants.R*products.solution.T))));

      for i in 1:nS loop
       uPure_substrates[i] = stateOfMatter.electroChemicalPotentialPure(
        substrates[i].definition,
        substrates[i].solution.T,
        substrates[i].solution.p,
        substrates[i].solution.v,
        substrates[i].solution.I);
      end for;

      for i in 1:nP loop
       uPure_products[i] = stateOfMatter.electroChemicalPotentialPure(
       products[i].definition,
       products[i].solution.T,
       products[i].solution.p,
       products[i].solution.v,
       products[i].solution.I);
      end for;

      Sx_fore = exp(s * ((substrates.state_forwards.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution.T)));
      Px_rear = exp((p * ((products.state_rearwards.u - uPure_products)./(Modelica.Constants.R*products.solution.T))));

      //debug
      Sx_rear = exp(s * ((substrates.state_rearwards.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution.T)));
      Px_fore = exp((p * ((products.state_forwards.u - uPure_products)./(Modelica.Constants.R*products.solution.T))));
      Kx = exp(- ((s * ((uPure_substrates)./(Modelica.Constants.R*substrates.solution.T))) - (p * ((uPure_products)./(Modelica.Constants.R*products.solution.T)))));

      //reaction molar rates
      rr*s = substrates.n_flow;
      rr*p = -products.n_flow;

      products.state_forwards.h = h_fore_mix*ones(nP);
      substrates.state_rearwards.h = h_rear_mix*ones(nS);

      if nS>0 then
        h_rear_mix*(substrates.n_flow*ones(nS)) + products.n_flow*products.state_rearwards.h = 0;
      else
        h_rear_mix = 0;
      end if;

      if nP>0 then
        h_fore_mix*(products.n_flow*ones(nP)) + substrates.n_flow*substrates.state_forwards.h = 0;
      else
        h_fore_mix = 0;
      end if;

      dr = (s * substrates.r) - (p * products.r);

      if nP>0 then

        (p * products.r) = (s * substrates.r)  -  der(rr)*L;

        for i in 2:nP loop
          //first product is based on inertial potential,
          //other products are provided as source with fixed flow and adaptation of their potential
          der(products[i].state_forwards.u).*TC = products[i].r;
        end for;
        for i in 2:nS loop
          der(substrates[i].state_rearwards.u).*TC = substrates[i].r;
        end for;
      end if;

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}})),
        Documentation(revisions="<html>
<p><i>2013-2025 by </i>Marek Mateják </p>
</html>",     info="<html>
<h4><span style=\"color: #008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>A<sub>i</sub></p></td>
<td><p>i-th substance</p></td>
</tr>
<tr>
<td><p>v<sub>i</sub></p></td>
<td><p>stochiometric coefficients of i-th substance</p></td>
</tr>
<tr>
<td><p>K</p></td>
<td><p>dissociation constant (activity based)</p></td>
</tr>
<tr>
<td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
<td><p>activity of the substance A</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub></p></td>
<td><p>activity coefficient of the substance A</p></td>
</tr>
<tr>
<td><p>x<sub>i</sub></p></td>
<td><p>mole fraction of the substance A</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
<td><p>molar enthalpy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
<td><p>molar Gibbs energy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
<td><p>molar entropy of formation of i-th substance</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>r</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by reaction</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
    end PartialReaction;
  end Internal;
  annotation (Documentation(info="<html>
<u>This package contains models implementing undirected versions of the processes. Here, the thermodynamic state of one or more fluid streams is changed by exchanging heat or work with the streams, or by delaying the state.</u>
</html>", revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>
</html>"), Icon(graphics={
         Ellipse(
          extent={{-60,54},{60,-66}},
          lineColor={158,66,200},
          lineThickness=0.5,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Line(
          points={{-94,0},{94,0}},
          color={158,66,200},
          thickness=0.5),
        Ellipse(
          extent={{-64,60},{56,-60}},
          lineColor={158,66,200},
          lineThickness=0.5,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end Processes;
