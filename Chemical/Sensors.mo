within Chemical;
package Sensors "Chemical sensors"
  extends Modelica.Icons.SensorsPackage;

  model MolarFlowSensor "Measure of molar flow"

    extends Modelica.Icons.RoundSensor;
    extends Interfaces.OnePort;

    Modelica.Blocks.Interfaces.RealOutput molarFlowRate(final unit="mol/s") annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-100})));

  equation
    molarFlowRate = port_a.q;

    port_a.u = port_b.u;

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

  model MoleFractionSensor "Measure of mole fraction"
    extends Modelica.Icons.RoundSensor;
    extends Interfaces.PartialSubstanceSensor;

    Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
    "Mole fraction of the substance"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={-100,0},
        rotation=180)));

  equation
    port_a.q = 0;

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

  model ElectroChemicalPotentialSensor
  "Measure of electro-chemical potential"
    extends Modelica.Icons.RoundSensor;

    Modelica.Blocks.Interfaces.RealOutput u(final unit="J/mol")
    "Electro-chemical potential of the substance"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={-100,0},
        rotation=180)));

  Interfaces.SubstancePort_b port_a
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  equation

    port_a.u = u;

    port_a.q = 0;
    port_a.h_outflow = 0;

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
  end ElectroChemicalPotentialSensor;

  model MolalitySensor "Measure of molality of the substance"
    extends Modelica.Icons.RoundSensor;
    extends Interfaces.PartialSubstanceSensor;

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
          origin={-100,0},
        rotation=180)));

  protected
  constant Modelica.Units.SI.Mass KG=1;
  equation
    port_a.q = 0;

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
    extends Interfaces.PartialSubstanceSensor;

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
          origin={-100,0},
        rotation=180)));

  protected
  constant Modelica.Units.SI.Volume L=0.001;
  equation
    port_a.q = 0;

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
    extends Interfaces.PartialSubstanceSensor;

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
          origin={-100,0},
        rotation=180)));

  equation
    port_a.q = 0;

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
    extends Interfaces.PartialSubstanceSensor;

     Modelica.Blocks.Interfaces.RealOutput partialPressure(final unit="Pa")
    "Partial pressure of the substance in gaseous solution"
     annotation (
        Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,-60}), iconTransformation(
          extent={{-20,-20},{20,20}},
          origin={-100,0},
        rotation=180)));

  equation
    port_a.q = 0;

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

    parameter Boolean useTemperatureInput = false
    "=true, if temperature is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Temperature T=system.T_ambient
    "Temperature if not useTemperatureInput" annotation (HideResult=
        true, Dialog(enable=not useTemperatureInput));

    Modelica.Blocks.Interfaces.RealInput temperature(start=
          T, final unit="K")=_temperature if useTemperatureInput
    "Temperature"
      annotation (HideResult=true,Placement(transformation(extent={{-120,58},
              {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40})));

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

  Interfaces.SubstancePort_b products[nP] "Products"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Interfaces.SubstancePort_b substrates[nS] "Substrates"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  Modelica.Units.SI.MolarEnergy DrG "Free Gibbs energy of reaction";

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
  Modelica.Units.SI.Temperature _temperature;
  Modelica.Units.SI.AmountOfSubstance _n;
  equation
    if not useTemperatureInput then
      _temperature = T;
    end if;
    if not useTotalAmountOfSubstancesInput then
      _n = n;
    end if;

    substrates.q = zeros(nS);
    substrates.h_outflow = zeros(nS);

    products.q = zeros(nP);
    products.h_outflow = zeros(nP);

    DrG = ((p * products.u) - (s * substrates.u));

    DissociationCoefficient_MoleFractionBased = exp(-DrG/(Modelica.Constants.R*T));

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

    parameter Boolean useTemperatureInput = false
    "=true, if temperature is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

  parameter Modelica.Units.SI.Temperature T=system.T_ambient
    "Temperature if not useTemperatureInput" annotation (HideResult=
        true, Dialog(enable=not useTemperatureInput));

    Modelica.Blocks.Interfaces.RealInput temperature(start=
          T, final unit="K")=_temperature if useTemperatureInput
    "Temperature"
      annotation (HideResult=true,Placement(transformation(extent={{-120,58},
              {-80,98}}), iconTransformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40})));

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

  Interfaces.SubstancePort_b products[nP] "Products"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Interfaces.SubstancePort_b substrates[nS] "Substrates"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  Modelica.Units.SI.MolarEnergy DrG "Free Gibbs energy of reaction";

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
  Modelica.Units.SI.Temperature _temperature;
  Modelica.Units.SI.AmountOfSubstance _n;
  equation
    if not useTemperatureInput then
      _temperature = T;
    end if;
    if not useTotalAmountOfSubstancesInput then
      _n = n;
    end if;

    substrates.q = zeros(nS);
    substrates.h_outflow = zeros(nS);

    products.q = zeros(nP);
    products.h_outflow = zeros(nP);

    DrG = ((p * products.u) - (s * substrates.u)) + (if (nP>0) then p[1] else 1)*Modelica.Constants.R*T*log(activityCoeficient);

    DissociationCoefficient_MoleFractionBased = exp(-DrG/(Modelica.Constants.R*T));

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
end Sensors;
