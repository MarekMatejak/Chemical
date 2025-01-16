within Chemical;
package Sensors "Chemical sensors"
  extends Modelica.Icons.SensorsPackage;

  model MolarFlowSensor "Measure of molar flow"

    extends Modelica.Icons.RoundSensor;

    Modelica.Blocks.Interfaces.RealOutput molarFlowRate(final unit="mol/s") annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-100})));

    Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
    Interfaces.Outlet outlet annotation (Placement(transformation(extent={{92,-10},{112,10}})));
  equation
    molarFlowRate = inlet.n_flow;

    connect(inlet, outlet) annotation (Line(
        points={{-98,0},{102,0}},
        color={158,66,200},
        thickness=0.5));
   annotation (
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Line(
            points={{70,-10},{90,-10}},
            color={127,0,127}),
          Line(
            points={{70,10},{90,10}},
            color={127,0,127}),
          Line(
            points={{-90,10},{-70,10}},
            color={127,0,127}),
          Line(
            points={{-90,-10},{-70,-10}},
            color={127,0,127}),
          Text(
            extent={{-31,-5},{28,-64}},
            lineColor={0,0,0},
            textString="dn")}));
  end MolarFlowSensor;

  model SubstanceMolarFlowSensor "Measure of molar flow with u0RT"

    extends Modelica.Icons.RoundSensor;

    Modelica.Blocks.Interfaces.RealOutput molarFlowRate(final unit="mol/s") annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-100})));

    Interfaces.InletProcess inlet annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
    Interfaces.OutletSubstance outlet annotation (Placement(transformation(extent={{92,-10},{112,10}})));
  equation
    molarFlowRate = inlet.n_flow;

    connect(inlet, outlet) annotation (Line(
        points={{-98,0},{102,0}},
        color={158,66,200},
        thickness=0.5));
   annotation (
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Line(
            points={{70,-10},{90,-10}},
            color={127,0,127}),
          Line(
            points={{70,10},{90,10}},
            color={127,0,127}),
          Line(
            points={{-90,10},{-70,10}},
            color={127,0,127}),
          Line(
            points={{-90,-10},{-70,-10}},
            color={127,0,127}),
          Text(
            extent={{-31,-5},{28,-64}},
            lineColor={0,0,0},
            textString="dn")}));
  end SubstanceMolarFlowSensor;

  model MoleFractionSensor "Measure of mole fraction"
    extends Modelica.Icons.RoundSensor;
    extends Internal.PartialSubstanceSensor;

    Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
    "Mole fraction of the substance"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  equation

    moleFraction = x;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
            textString="x"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end MoleFractionSensor;

  model SubstanceMoleFractionSensor "Measure of mole fraction"
    extends Modelica.Icons.RoundSensor;

   Interfaces.InletProcess inlet "The substance" annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

    Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
    "Mole fraction of the substance"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  equation

    inlet.n_flow = 0;

    moleFraction = exp(inlet.uRT - inlet.u0RT);

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
            textString="x"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end SubstanceMoleFractionSensor;

  model URTSensor
    "Measure of electro-chemical potential divided by gas constant and temperature"
    extends Modelica.Icons.RoundSensor;

    Modelica.Blocks.Interfaces.RealOutput uRT(final unit="J/mol")
    "Electro-chemical potential of the substance divided by gas constant and temperature"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={-100,0},
        rotation=180)));

    Interfaces.Inlet port_a annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  equation

    port_a.uRT = uRT;

    port_a.n_flow = 0;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
          textString="u"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end URTSensor;

  model MolalitySensor "Measure of molality of the substance"
    extends Modelica.Icons.RoundSensor;
    extends Internal.PartialSubstanceSensor;

  parameter Modelica.Units.SI.AmountOfSubstance
    AmountOfSolutionPer1kgOfSolvent=1
    "Amount of all particles in the solution per one kilogram of solvent";

     Modelica.Blocks.Interfaces.RealOutput molality(final unit="mol/kg")
    "Molality of the substance (amount of substance per mass of solvent)"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  protected
  constant Modelica.Units.SI.Mass KG=1;
  equation

    x=molality*KG / AmountOfSolutionPer1kgOfSolvent;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
          textString="b"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end MolalitySensor;

  model MolarConcentrationSensor "Measure of molarity of the substance"
    extends Modelica.Icons.RoundSensor;
    extends Internal.PartialSubstanceSensor;

  parameter Modelica.Units.SI.AmountOfSubstance
    AmountOfSolutionInOneLiter=1
    "Amount of all particles in one liter of the solution";

     Modelica.Blocks.Interfaces.RealOutput molarConcentration(final unit="mol/m3", displayUnit="mol/l")
    "Molarity of the substance (amount of substance in one liter of whole solution)"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  protected
  constant Modelica.Units.SI.Volume L=0.001;
  equation

    x=molarConcentration*L / AmountOfSolutionInOneLiter;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
          textString="c"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end MolarConcentrationSensor;

  model MassFractionSensor "Measure of mass fraction of the substance"
    extends Modelica.Icons.RoundSensor;
    extends Internal.PartialSubstanceSensor;

  parameter Modelica.Units.SI.AmountOfSubstance
    AmountOfSolutionInOneKilogram=1
    "Amount of all particles in one kilogram of the solution";

     Modelica.Blocks.Interfaces.RealOutput massFraction(final unit="kg/kg")
    "Mass fraction of the substance (mass of the substance per mass of the whole solution)"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  equation

    x=(massFraction*stateOfMatter.specificAmountOfParticles(substanceData)) / AmountOfSolutionInOneKilogram;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
          textString="mx"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end MassFractionSensor;

  model PartialPressureSensor
  "Measure of partial pressure of the substance in gaseous solution"
    extends Modelica.Icons.RoundSensor;
    extends Internal.PartialSubstanceSensor;

     Modelica.Blocks.Interfaces.RealOutput partialPressure(final unit="Pa")
    "Partial pressure of the substance in gaseous solution"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={100,0},
        rotation=0)));

  equation

    partialPressure = x*solution.p;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-31,-3},{28,-62}},
            lineColor={0,0,0},
          textString="p"),
          Line(
            points={{70,0},{80,0}},
            color={127,0,127})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end PartialPressureSensor;

  model DissociationCoefficient
  "Meassure dissociation coefficient (mole fraction based) for pure substances"
    extends Modelica.Icons.RectangularSensor;

    outer Modelica.Fluid.System system "System wide properties";



    parameter Boolean useTotalAmountOfSubstancesInput = false
    "=true, if total amount of substances in solution is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

  parameter Modelica.Units.SI.AmountOfSubstance n=1
    "Amount of all substances in solution per one liter of solution if not useTotalAmountOfSubstancesInput"
    annotation (HideResult=true, Dialog(enable=not
          useTotalAmountOfSubstancesInput));

    Modelica.Blocks.Interfaces.RealInput totalAmountOfSubstances(start=
          n, final unit="mol")=_n if useTotalAmountOfSubstancesInput
    "Temperature"
      annotation (HideResult=true,Placement(transformation(extent={{-120,58},
              {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={40,40})));

  parameter Modelica.Units.SI.Mass m=1
    "Mass of solvent per one liter of solution";

    parameter Integer nS=0 "Number of substrates types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

  parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
    "Stoichiometric reaction coefficient for substrates"
    annotation (HideResult=true);

    parameter Integer nP=0 "Number of products types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

  parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
    "Stoichiometric reaction coefficients for products"
    annotation (HideResult=true);

    Interfaces.Outlet                   products[nP] "Products" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Interfaces.Inlet                    substrates[nS] "Substrates" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    Chemical.Utilities.Units.URT GRT "Free Gibbs energy of reaction divided by gas constant and temperature";

    Modelica.Blocks.Interfaces.RealOutput DissociationCoefficient_MoleFractionBased
    "Dissociation constant (if all substances has activity=1)"   annotation (Placement(transformation(
            extent={{-6,-86},{14,-66}}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-80})));

    Real DissociationCoefficient_MolalityBased
    "As ratio of molalities in moles per 1 kg of solvent";
    Real DissociationCoefficient_MolarityBased
    "As ratio of molar concentration in moles per liter of solution";

    Real pK
    "= -log10('mole-fraction based dissociation coefficient')";

  protected
  Modelica.Units.SI.AmountOfSubstance _n;
  equation

    if not useTotalAmountOfSubstancesInput then
      _n = n;
    end if;

    substrates.n_flow = zeros(nS);


    products.n_flow = zeros(nP);
    products.h = zeros(nP);

    GRT = ((p * products.uRT) - (s * substrates.uRT));

    DissociationCoefficient_MoleFractionBased = exp(-GRT);

    pK=-log10(DissociationCoefficient_MoleFractionBased);

    DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

    DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-160,-94},{-12,-68}},
            lineColor={0,0,0},
          textString="%s"),
          Text(
            extent={{12,-92},{160,-66}},
            lineColor={0,0,0},
          textString="%p")}),
      Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
  end DissociationCoefficient;

  model ActivityCoefficient
    "Calculate activity coefficient for product[1]"
    extends Modelica.Icons.RectangularSensor;

    outer Modelica.Fluid.System system "System wide properties";

    parameter Boolean useTotalAmountOfSubstancesInput = false
    "=true, if total amount of substances in solution is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

  parameter Modelica.Units.SI.AmountOfSubstance n=1
    "Amount of all substances in solution per one liter of solution if not useTotalAmountOfSubstancesInput"
    annotation (HideResult=true, Dialog(enable=not
          useTotalAmountOfSubstancesInput));

    Modelica.Blocks.Interfaces.RealInput totalAmountOfSubstances(start=
          n, final unit="mol")=_n if useTotalAmountOfSubstancesInput
    "Temperature"
      annotation (HideResult=true,Placement(transformation(extent={{-120,58},
              {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={40,40})));

  parameter Modelica.Units.SI.Mass m=1
    "Mass of solvent per one liter of solution";

    parameter Integer nS=0 "Number of substrates types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

  parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
    "Stoichiometric reaction coefficient for substrates"
    annotation (HideResult=true);

    parameter Integer nP=0 "Number of products types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, group="Ports"));

  parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
    "Stoichiometric reaction coefficients for products"
    annotation (HideResult=true);

    Interfaces.Outlet                   products[nP] "Products" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Interfaces.Inlet                    substrates[nS] "Substrates" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  Chemical.Utilities.Units.URT GRT "Free Gibbs energy of reaction divided by gas constant and temperature";


    Modelica.Blocks.Interfaces.RealOutput activityCoeficient
    "Activity coeficient of one product"   annotation (Placement(transformation(
            extent={{-6,-86},{14,-66}}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-80})));

    parameter Boolean MolarityBased = true "if dissociation coefficient is molarity based";

    parameter Real DissociationCoefficient_MoleFractionBased = if MolarityBased then DissociationCoefficient_MolarityBased/((n/1)^(p*ones(nP)-s*ones(nS))) else DissociationCoefficient_MolalityBased/((n/m)^(p*ones(nP)-s*ones(nS)))
    "K as ratio of mole fractions";
    parameter Real DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased
    "K as ratio of molalities in moles per 1 kg of solvent"
    annotation (HideResult=true, Dialog(enable=not MolarityBased));
    parameter Real DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased
    "K as ratio of molar concentration in moles per liter of solution"
    annotation (HideResult=true, Dialog(enable=MolarityBased));

    Real pK
    "= -log10('mole-fraction based dissociation coefficient')";

  protected
  Modelica.Units.SI.AmountOfSubstance _n;
  equation
    if not useTotalAmountOfSubstancesInput then
      _n = n;
    end if;

    substrates.n_flow = zeros(nS);

    products.n_flow = zeros(nP);
    products.h = zeros(nP);

    GRT = ((p * products.uRT) - (s * substrates.uRT)) + (if (nP>0) then p[1] else 1)*log(activityCoeficient);

    DissociationCoefficient_MoleFractionBased = exp(-GRT);

    pK=-log10(DissociationCoefficient_MoleFractionBased);

    //DissociationCoefficient_MolalityBased = ((n/m)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

    //DissociationCoefficient_MolarityBased = ((n/1)^(p*ones(nP)-s*ones(nS))) * DissociationCoefficient_MoleFractionBased;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Text(
            extent={{-160,-94},{-12,-68}},
            lineColor={0,0,0},
          textString="%s"),
          Text(
            extent={{12,-92},{160,-66}},
            lineColor={0,0,0},
          textString="%p")}),
      Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
  end ActivityCoefficient;

  package Internal
    partial model PartialSubstance

     outer Modelica.Fluid.System system "System wide properties";

      Interfaces.Inlet inlet "The substance" annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

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

     parameter stateOfMatter.SubstanceData substanceData
     "Definition of the substance"
        annotation (choicesAllMatching = true);

    Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

    Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

    protected
    Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

    Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";

    Modelica.Units.SI.Temperature temperature
      "Temperature of the solution";

    Modelica.Units.SI.Pressure pressure "Pressure of the solution";

    Modelica.Units.SI.ElectricPotential electricPotential
      "Electric potential of the solution";

    Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength
      "Ionic strength of the solution";

    //Modelica.Units.SI.MolarMass molarMass "Molar mass of the substance";

    Modelica.Units.SI.MolarEnthalpy molarEnthalpy
      "Molar enthalpy of the substance";

    Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";

    Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

    Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

    Modelica.Units.SI.MolarVolume molarVolume
      "Molar volume of the substance";

    Modelica.Units.SI.MolarVolume molarVolumePure
      "Molar volume of the pure substance";

    Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

      //  Modelica.SIunits.MolarHeatCapacity molarHeatCapacityCp
      //    "Molar heat capacity of the substance at constant pressure";

    equation
     //aliases
     gamma = stateOfMatter.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     z = stateOfMatter.chargeNumberOfIon(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
    // molarMass = stateOfMatter.molarMass(substanceData);

     molarEnthalpy = stateOfMatter.molarEnthalpy(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarEntropyPure = stateOfMatter.molarEntropyPure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     u0 = stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
     uPure = stateOfMatter.electroChemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
     molarVolume = stateOfMatter.molarVolume(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarVolumePure = stateOfMatter.molarVolumePure(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     molarVolumeExcess = stateOfMatter.molarVolumeExcess(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
     //  molarHeatCapacityCp = stateOfMatter.molarHeatCapacityCp(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     inlet.uRT = stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength)/(Modelica.Constants.R*temperature)
       + log(a)
       + z*Modelica.Constants.F*electricPotential/(Modelica.Constants.R*temperature);

     inlet.h = molarEnthalpy;

     annotation (
       Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialSubstance;

    partial model PartialSubstanceInSolution "Substance properties for components, where the substance is connected with the solution"

      Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

      extends PartialSubstance;

    protected
    Modelica.Units.SI.AmountOfSubstance amountOfSolution
      "Amount of all solution particles";

    equation

      //aliases
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;

    end PartialSubstanceInSolution;

    model PartialSubstanceSensor "Base class for sensor based on substance and solution properties"
      extends PartialSubstanceInSolution;

    equation
      //solution is not changed by the sensor components
      solution.dH = 0;
      solution.i = 0;
      solution.dV = 0;
      solution.Gj = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.Vj = 0;

    end PartialSubstanceSensor;

    package Types
      type Quantities = enumeration(
          c_molpm3 "Concentration (mmol/L)",
          X_kgpkg "Mass fraction (kg/kg)",
          b_molpkg "Molality (mol/kg)",
          x_molpmol "Mole fraction (mol/mol)",
          p_Pa "Partial pressure (Pa)",
          p_mmHg "Partial pressure (mmHg)",
          p_bar "Partial pressure (bar)",
          u_Jpmol "Steadystate electro-chemical potential (J/mol)",
          u_kJpmol "Steadystate electro-chemical potential (kJ/mol)",
          r_Jpmol "Inertial electro-chemical potential (J/mol)",
          r_kJpmol "Inertial electro-chemical potential (kJ/mol)",
          u_total_Jpmol "Total electro-chemical potential (J/mol)",
          u_total_kJpmol "Total electro-chemical potential (kJ/mol)",
          h_Jpmol "Specific enthalpy (J/mol)",
          s_JpmolK "Specific enthropy (J/(mol.K))");
      type InitializationModelSensor = enumeration(
        steadyState
          "Steady state initialization (derivatives of states are zero)",
        state
          "Initialization with initial output state") "Initialization modes for sensor lowpass";
      function getQuantity "Computes selected quantity from state"
        extends Modelica.Icons.Function;

        replaceable package Medium =
            ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
          "Medium model"
          annotation (choicesAllMatching=true,
            Documentation(info="<html>
      <p>Medium Model for the function. Make sure it implements the needed functions.</p>
        </html>"));

        input Medium.ThermodynamicState state;
        input Modelica.Units.SI.Pressure r;
        input ThermofluidStream.Sensors.Internal.Types.Quantities quantity;
        input Modelica.Units.SI.Density rho_min;
        output Real value;

      algorithm
        if quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.T_K then
          value := Medium.temperature(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.T_C then
          value :=Modelica.Units.Conversions.to_degC(Medium.temperature(state));
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.p_Pa then
          value := Medium.pressure(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar then
          value :=Modelica.Units.Conversions.to_bar(Medium.pressure(state));
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.r_Pa then
          value := r;
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.r_bar then
          value :=Modelica.Units.Conversions.to_bar(r);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.p_total_Pa then
          value := Medium.pressure(state)+r;
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.p_total_bar then
          value :=Modelica.Units.Conversions.to_bar(Medium.pressure(state) + r);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.h_Jpkg then
          value := Medium.specificEnthalpy(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.s_JpkgK then
          value := Medium.specificEntropy(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.rho_kgpm3 then
          value := Medium.density(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.v_m3pkg then
          value := 1/(max(rho_min, Medium.density(state)));
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.a_mps then
          value := Medium.velocityOfSound(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.cv_JpkgK then
          value := Medium.specificHeatCapacityCv(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.cp_JpkgK then
          value := Medium.specificHeatCapacityCp(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.kappa_1 then
          value := Medium.isentropicExponent(state);
        elseif quantity == ThermofluidStream.Sensors.Internal.Types.Quantities.MM_kgpmol then
          value := Medium.molarMass(state);
        else
          value :=0;
        end if;

        annotation (Documentation(info="<html>
<p>Helper function to get a quantity from an Thermofluid state.</p>
</html>"));
      end getQuantity;
    end Types;

    function getQuantity "Computes selected quantity from state"
      extends Modelica.Icons.Function;

      replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
        Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
        annotation (choicesAllMatching=true,
          Documentation(info="<html>
      <p>Medium Model for the function. Make sure it implements the needed functions.</p>
        </html>"));

      input Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential";
      input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";
      input Modelica.Units.SI.Pressure r "Inertial electro-chemical potential";
      input Types.Quantities quantity "What to measure?";
      input stateOfMatter.SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature temperature=298.15 "Temperature";
      input Modelica.Units.SI.Pressure pressure=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential electricPotential=0
       "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0
       "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass solutionMass;
      input Modelica.Units.SI.AmountOfSubstance solutionAmount;
      input Modelica.Units.SI.Volume solutionVolume;

      output Real value;

      /*
    c_molpm3 "Concentration (mmol/L)",
    X_kgpkg "Mass fraction (kg/kg)",
    b_molpkg "Molality (mol/kg)",
    x_molpmol "Mole fraction (mol/mol)",
    p_Pa "Partial pressure (Pa)",
    p_mmHg "Partial pressure (mmHg)",
    p_bar "Partial pressure (bar)",
    u_Jpmol "Steadystate electro-chemical potential (J/mol)",
    u_kJpmol "Steadystate electro-chemical potential (kJ/mol)",
    r_Jpmol "Inertial electro-chemical potential (J/mol)",
    r_kJpmol "Inertial electro-chemical potential (kJ/mol)",
    u_total_Jpmol "Total electro-chemical potential (J/mol)",
    u_total_kJpmol "Total electro-chemical potential (kJ/mol)",
    h_Jpmol "Specific enthalpy (J/mol)",
    s_JpmolK "Specific enthropy (J/(mol.K))"
    */

    protected
      Modelica.Units.SI.ChargeNumberOfIon z;
      Modelica.Units.SI.ChemicalPotential u0;
      Modelica.Units.SI.MoleFraction a,x;
      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

    algorithm
      z :=  stateOfMatter.chargeNumberOfIon( substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
      gamma := stateOfMatter.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);


      u0 := stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength)
       + z*Modelica.Constants.F*electricPotential;

      a := exp((u - u0)/(Modelica.Constants.R*temperature));
      x := a/gamma;

      if quantity == Types.Quantities.c_molpm3 then
        value := (x * solutionAmount)/solutionVolume;
      elseif quantity == Types.Quantities.X_kgpkg then
        value := ((x * solutionAmount)/solutionMass)/stateOfMatter.specificAmountOfParticles(substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength);
      elseif quantity == Types.Quantities.b_molpkg then
        value := (x * solutionAmount)/solutionMass;
      elseif quantity == Types.Quantities.x_molpmol then
        value := x;
      elseif quantity == Types.Quantities.p_Pa then
        value := x*pressure;
      elseif quantity == Types.Quantities.p_mmHg then
        value := x*pressure * (760/101325);
      elseif quantity == Types.Quantities.p_bar then
        value := Modelica.Units.Conversions.to_bar(x*pressure);
      elseif quantity == Types.Quantities.u_Jpmol then
        value := u;
      elseif quantity == Types.Quantities.u_kJpmol then
        value := u/1000;
      elseif quantity == Types.Quantities.h_Jpmol then
        value := h;
      elseif quantity == Types.Quantities.s_JpmolK then
        value := (h-u)/temperature;
      else
        value :=0;
      end if;

      annotation (Documentation(info="<html>
<p>Helper function to get a quantity from an Thermofluid state.</p>
</html>"));
    end getQuantity;

    function getUnit "Returns unit of input quantity"
      extends Modelica.Icons.Function;

      input Types.Quantities quantity;
      output String unit;

    algorithm

      if quantity == Types.Quantities.c_molpm3 then
        unit := "mol/m3";
      elseif quantity == Types.Quantities.X_kgpkg then
        unit := "kg/kg";
      elseif quantity == Types.Quantities.b_molpkg then
        unit := "mol/kg";
      elseif quantity == Types.Quantities.x_molpmol then
        unit := "mol/mol";
      elseif quantity == Types.Quantities.p_Pa then
        unit := "Pa";
      elseif quantity == Types.Quantities.p_mmHg then
        unit := "mmHg";
      elseif quantity == Types.Quantities.p_bar then
        unit := "bar";
      elseif quantity == Types.Quantities.u_Jpmol then
        unit := "J/mol";
      elseif quantity == Types.Quantities.u_kJpmol then
        unit := "kJ/mol";
      elseif quantity == Types.Quantities.h_Jpmol then
        unit :="J/mol";
      elseif quantity == Types.Quantities.s_JpmolK then
        unit := "J/(mol.K)";
      else
        unit :="";
      end if;

      annotation (Documentation(info="<html>
<p>Helper function to get the unit for a quantity.</p>
</html>"));
    end getUnit;
  end Internal;

  model SingleSensorSelect "Sensor with selectable measured quantity"
    import Chemical.Sensors.Internal.Types.Quantities;
    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;
    extends Internal.PartialSubstanceSensor;

    parameter Integer digits(min=0) = 1 "Number of displayed digits";
    parameter Quantities quantity "Quantity the sensor measures";
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real value_0(unit=Internal.getUnit(quantity)) = 0 "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init == InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

    Modelica.Blocks.Interfaces.RealOutput value_out(unit=Internal.getUnit(quantity)) = value if outputValue "Measured value [variable]"
      annotation (Placement(transformation(extent={{80,-20},{120,20}})));

    output Real value(unit=Internal.getUnit(quantity)) "Computed value of the selected quantity";

  protected
    outer DropOfCommons dropOfCommons;

    Real direct_value(unit=Internal.getUnit(quantity));

    function getQuantity = Internal.getQuantity (redeclare package
          stateOfMatter =
            stateOfMatter)                                                       "Quantity compute function"
      annotation (Documentation(info="<html>
      <p>This function computes the selected quantity from state. r and rho_min are neddet for the quantities r/p_total and v respectively.</p>
      </html>"));

  initial equation
    if filter_output and init==InitMode.steadyState then
      value= direct_value;
    elseif filter_output then
      value = value_0;
    end if;

  equation

    direct_value = getQuantity(inlet.u, inlet.h, inlet.r, quantity, substanceData, temperature, pressure, electricPotential, moleFractionBasedIonicStrength,
     solution.m, solution.n, solution.V);


    if filter_output then
      der(value) * TC = direct_value-value;
    else
      value = direct_value;
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-54,24},{66,-36}},
            lineColor={0,0,0},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-100,0},{0,0}},
            color={28,108,200},
            thickness=0.5),
          Rectangle(
            extent={{-60,30},{60,-30}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-60,30},{60,-30}},
            textColor={28,108,200},
            textString=DynamicSelect("value", String(
                value,
                format="1."+String(digits)+"f"))),
          Text(
            extent={{0,25},{60,75}},
            textColor={175,175,175},
            textString="%quantity")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Sensor for measuring a selectable quantity.</p>
<p>This sensor can be connected to a fluid stream without a junction.</p>
</html>"));
  end SingleSensorSelect;
end Sensors;
