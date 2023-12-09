within Chemical;
package Sources "Chemical sources"
  extends Modelica.Icons.SourcesPackage;

  model PureSubstance "Constant source of pure substance"
    extends Interfaces.PartialSubstance;

    parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
      "Temperature"
      annotation (HideResult=true, Dialog(enable=not useTemperatureInput));
    parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
      "Pressure"
      annotation (HideResult=true, Dialog(enable=not usePressureInput));
    parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
      "Electric potential" annotation (HideResult=true, Dialog(enable=not
            useElectricPotentialInput));
    parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
       0 "Ionic strength" annotation (HideResult=true, Dialog(enable=not
            useIonicStrengthInput));

    parameter Boolean useTemperatureInput = false
    "=true, if temperature is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Boolean usePressureInput = false
    "=true, if pressure is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Boolean useElectricPotentialInput = false
    "=true, if electric potential is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Boolean useIonicStrengthInput = false
    "=true, if ionic strength is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput T(start=
          Temperature, final unit="K")=temperature if useTemperatureInput
    "Temperature"
      annotation (HideResult=true,Placement(transformation(extent={{-120,58},
              {-80,98}}), iconTransformation(extent={{-120,58},{-80,98}})));

    Modelica.Blocks.Interfaces.RealInput p(start=
          Pressure, final unit="Pa")=pressure if usePressureInput
    "Pressure"
      annotation (HideResult=true,Placement(transformation(extent={{-120,16},
              {-80,56}}), iconTransformation(extent={{-120,16},{-80,56}})));

    Modelica.Blocks.Interfaces.RealInput v(start=
          ElectricPotential, final unit="Pa")=electricPotential if useElectricPotentialInput
    "Electric potential"
      annotation (HideResult=true,Placement(transformation(extent={{-120,-60},
              {-80,-20}}),iconTransformation(extent={{-120,-60},{-80,-20}})));

    Modelica.Blocks.Interfaces.RealInput I(start=
          MoleFractionBasedIonicStrength, final unit="mol/mol")=moleFractionBasedIonicStrength if useIonicStrengthInput
    "Pressure"
      annotation (HideResult=true,Placement(transformation(extent={{-120,-100},
              {-80,-60}}),iconTransformation(extent={{-120,-100},{-80,-60}})));
  protected
    Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
        Modelica.Constants.R*temperature))
      "Dissociation constant of hydrogen bond between base molecules";
    Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
        stateOfMatter.selfClusteringBondEnthalpy(
                                             substanceData) - temperature*
        stateOfMatter.selfClusteringBondEntropy(
                                            substanceData)
      "Gibbs energy of hydrogen bond between H2O molecules";

  equation

     if stateOfMatter.selfClustering(substanceData) then

      //Liquid cluster theory - equilibrium:
      //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules

      //sum(x[i]) = x/(1-K*x) = amountOfParticles/amountOfParticles = 1;
      x = 1/(1+SelfClustering_K) "mole fraction of free base molecule";
    else
      x = 1 "pure substance is composed only with free base molecules";
    end if;

     if (not useTemperatureInput) then
       temperature = Temperature;
     end if;
     if (not usePressureInput) then
       pressure = Pressure;
     end if;
     if (not useElectricPotentialInput) then
       electricPotential = ElectricPotential;
     end if;
     if (not useIonicStrengthInput) then
       moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;
     end if;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={107,45,134},
            fillPattern=FillPattern.Backward),
          Text(
            extent={{10,8},{-90,-92}},
            lineColor={0,0,0},
            textString="pure"),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-104,-76},{100,-100}},
            lineColor={0,0,0},
            textString="%T K")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end PureSubstance;

  model ExternalIdealGasSubstance
  "Ideal gas substance with defined partial pressure"
    extends Interfaces.PartialSubstance(redeclare package stateOfMatter =
          Interfaces.IdealGas);

    parameter Boolean usePartialPressureInput = false
    "=true, if fixed partial pressure is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Pressure PartialPressure=0
      "Fixed partial pressure if usePartialPressureInput=false" annotation (
       HideResult=true, Dialog(enable=not usePartialPressureInput));

    parameter Modelica.Units.SI.Pressure TotalPressure=system.p_ambient
      "Total pressure of the whole gaseous solution";

    parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
      "Temperature";
    parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
       0 "Ionic strength";
    parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
      "Electric potential";

    Modelica.Blocks.Interfaces.RealInput partialPressure(start=
          PartialPressure, final unit="Pa")=p if usePartialPressureInput
    "Partial pressure of gas = total pressure * gas fraction"
      annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.Pressure p "Current partial pressure";

    parameter Modelica.Units.SI.Volume Volume=0.001
      "Volume of gaseous solution";

  equation
    if not usePartialPressureInput then
      p=PartialPressure;
    end if;

    //mole fraction
    x = p / TotalPressure;

    //the solution
    temperature = Temperature;
    pressure = TotalPressure;
    electricPotential = ElectricPotential;
    moleFractionBasedIonicStrength = 0;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          pattern=LinePattern.None,
          fillColor={170,255,255},
          fillPattern=FillPattern.Backward),
          Polygon(
            points={{-100,100},{100,-100},{100,100},{-100,100}},
            fillColor={159,159,223},
            fillPattern=FillPattern.Backward,
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Text(
            extent={{0,0},{-100,-100}},
            lineColor={0,0,0},
            textString="P,T"),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-100,-102},{104,-126}},
            lineColor={0,0,0},
            textString="%T K")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ExternalIdealGasSubstance;

  model ExternalMolality "Constant source of substance molality"
    extends Interfaces.PartialSubstance;

    outer Modelica.Fluid.System system "System wide properties";

     parameter Real Molality(final unit="mol/kg") = 1e-8
    "Fixed molality of the substance if useMolalityInput=false"
      annotation (HideResult=true, Dialog(enable=not useMolalityInput));

    parameter Modelica.Units.SI.AmountOfSubstance
      AmountOfSolutionPer1KgSolvent=55.508
      "Amount of all particles in the solution per one kilogram of solvent";

      parameter Boolean useMolalityInput = false
    "Is amount of substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
      "Temperature";
    parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
      "Pressure";
    parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
       0 "Ionic strength";
    parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
      "Electric potential";

    Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG
      if useMolalityInput
      annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

  protected
    constant Modelica.Units.SI.Mass KG=1;
  equation
     if not useMolalityInput then
       n=Molality*KG;
     end if;

    x = n/AmountOfSolutionPer1KgSolvent;

    //solution properties at the port
    temperature = Temperature;
    pressure = Pressure;
    electricPotential = ElectricPotential;
    moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={107,45,134},
            fillPattern=FillPattern.Backward),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-104,-76},{100,-100}},
            lineColor={0,0,0},
            textString="%T K"),
          Text(
            extent={{94,-4},{-94,-78}},
            lineColor={0,0,0},
            textString="molality")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ExternalMolality;

  model ExternalConcentration "Constant source of molar concentration"
     extends Interfaces.PartialSubstance;

     outer Modelica.Fluid.System system "System wide properties";

     parameter Real MolarConcentration(final unit="mol/m3", displayUnit="mol/l") = 1e-8
    "Fixed molarity of the substance if useMolarityInput=false"
      annotation (HideResult=true, Dialog(enable=not useMolarityInput));

    parameter Modelica.Units.SI.AmountOfSubstance AmountOfSolutionIn1L=55.508
      "Amount of all particles in the solution one liter of solvent";

      parameter Boolean useMolarityInput = false
    "Is amount of substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
      "Temperature";
    parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
      "Pressure";
    parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
       0 "Ionic strength";
    parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
      "Electric potential";

    Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L
      if useMolarityInput
      annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

  protected
    constant Modelica.Units.SI.Volume L=0.001;
  equation
     if not useMolarityInput then
       n=MolarConcentration*L;
     end if;

    x = n/AmountOfSolutionIn1L;

    //solution properties at the port
    temperature = Temperature;
    pressure = Pressure;
    electricPotential = ElectricPotential;
    moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={107,45,134},
            fillPattern=FillPattern.Backward),
          Text(
            extent={{94,92},{-94,18}},
            lineColor={0,0,0},
            textString="molarity"),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-104,-76},{100,-100}},
            lineColor={0,0,0},
            textString="%T K")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ExternalConcentration;

  model ExternalMoleFraction "Constant source of substance mole fraction"
       extends Interfaces.PartialSubstance;

     outer Modelica.Fluid.System system "System wide properties";

    parameter Modelica.Units.SI.MoleFraction MoleFraction=1e-8
      "Fixed mole fraction of the substance if useMoleFractionInput=false"
      annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

      parameter Boolean useMoleFractionInput = false
    "Is mole fraction of the substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Temperature Temperature=system.T_ambient
      "Temperature";
    parameter Modelica.Units.SI.Pressure Pressure=system.p_ambient
      "Pressure";
    parameter Modelica.Units.SI.MoleFraction MoleFractionBasedIonicStrength=
       0 "Ionic strength";
    parameter Modelica.Units.SI.ElectricPotential ElectricPotential=0
      "Electric potential";

    Modelica.Blocks.Interfaces.RealInput moleFractionInput(
      final unit="mol/mol",
      start=MoleFraction)=x
      if useMoleFractionInput annotation (HideResult=true, Placement(transformation(
            extent={{-120,-20},{-80,20}})));

  equation
     if not useMoleFractionInput then
       x=MoleFraction;
     end if;

    //solution properties at the port
    temperature = Temperature;
    pressure = Pressure;
    electricPotential = ElectricPotential;
    moleFractionBasedIonicStrength = MoleFractionBasedIonicStrength;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={107,45,134},
            fillPattern=FillPattern.Backward),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-104,-76},{100,-100}},
            lineColor={0,0,0},
            textString="%T K"),
          Text(
            extent={{94,-4},{-94,-78}},
            lineColor={0,0,0},
            textString="n")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ExternalMoleFraction;

  model ExternalElectroChemicalPotential
  "Constant source of electro-chemical potential"

  parameter Modelica.Units.SI.ChemicalPotential U=1e-8
    "Fixed electro-chemical potential of the substance if usePotentialInput=false"
    annotation (HideResult=true, Dialog(enable=not usePotentialInput));

     parameter Boolean usePotentialInput = false
    "Is electro-chemical potential of the substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput uInput(final unit="J/mol")=port_a.u
      if usePotentialInput annotation (HideResult=true, Placement(transformation(
            extent={{-120,-20},{-80,20}})));

  Interfaces.SubstancePort_a port_a
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
  equation
     if not usePotentialInput then
       port_a.u=U;
     end if;

    port_a.h_outflow = MolarHeat;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={107,45,134},
            fillPattern=FillPattern.Backward),
          Line(
            points={{-62,0},{56,0}},
            color={191,0,0},
            thickness=0.5),
          Polygon(
            points={{38,-20},{38,20},{78,0},{38,-20}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-104,-76},{100,-100}},
            lineColor={0,0,0},
            textString="%T K"),
          Text(
            extent={{94,-4},{-94,-78}},
            lineColor={0,0,0},
            textString="molality")}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ExternalElectroChemicalPotential;

  model SubstanceInflow "Molar pump of substance to system"
    extends Interfaces.ConditionalSubstanceFlow;

  Interfaces.SubstancePort_b port_b "Outflow"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
  equation
    port_b.q = -q;
    port_b.h_outflow = MolarHeat;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics={
          Rectangle(
            extent={{-100,-42},{100,40}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-48,20},{50,0},{-48,-21},{-48,20}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-220,-80},{220,-60}},
            textString="%name",
            lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end SubstanceInflow;

  model SubstanceOutflow "Molar pump of substance out of system"
    extends Interfaces.ConditionalSubstanceFlow;

  Interfaces.SubstancePort_b port_a "Inflow"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  parameter Modelica.Units.SI.MolarEnthalpy MolarHeat=0;
  equation
    port_a.q = q;

    port_a.h_outflow = MolarHeat;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics={
          Rectangle(
            extent={{-100,-42},{100,40}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-48,20},{50,0},{-48,-21},{-48,20}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-220,-80},{220,-60}},
            textString="%name",
            lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end SubstanceOutflow;

  model Clearance "Physiological Clearance"
   extends Boundaries.Internal.ConditionalSolutionFlow(
                                              final SolutionFlow=Clearance/K);
   extends Interfaces.PartialSubstanceSensor;

  parameter Modelica.Units.SI.VolumeFlowRate Clearance=0
    "Physiological clearance of the substance if useSolutionFlowInput=false"
    annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

    parameter Real K(unit="1")=1
    "Coefficient such that Clearance = K*solutionFlow";

  Modelica.Units.SI.MolarFlowRate molarClearance
    "Current molar clearance";

  equation
    molarClearance = q*K;

    port_a.q = molarClearance * x;

    assert(molarClearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                          graphics={
          Rectangle(
            extent={{-100,-100},{100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{80,25},{-80,0},{80,-25},{80,25}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-90},{150,-50}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-100,-30},{100,-50}},
            lineColor={0,0,0},
            textString="K=%K")}),        Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Clearance;

  model Degradation "Degradation of the substance"
    extends Interfaces.PartialSubstanceSensor;

  parameter Modelica.Units.SI.Time HalfTime
    "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

  equation
    port_a.q = (Modelica.Math.log(2)/HalfTime)*x*amountOfSolution;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                          graphics={
          Rectangle(
            extent={{-100,-100},{100,58}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{64,26},{-78,0},{64,-26},{64,26}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-148,-82},{152,-42}},
            textString="%name",
            lineColor={128,0,255}),
          Text(
            extent={{-100,54},{100,28}},
            lineColor={0,0,0},
            textString="t1/2 = %HalfTime s"),
          Polygon(
            points={{54,24},{54,-24},{44,-22},{44,22},{54,24}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{30,20},{30,-20},{20,-18},{20,18},{30,20}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{8,16},{8,-16},{-2,-14},{-2,14},{8,16}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,12},{-12,-12},{-22,-10},{-22,10},{-12,12}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-34,8},{-34,-8},{-44,-6},{-44,6},{-34,8}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-56,4},{-56,-4},{-66,-2},{-66,2},{-56,4}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)}),
      Documentation(revisions="<html>
<table>
<tr>
<td>Author:</td>
<td>Marek Matejak</td>
</tr>
<tr>
<td>Copyright:</td>
<td>In public domains</td>
</tr>
<tr>
<td>By:</td>
<td>Charles University, Prague</td>
</tr>
<tr>
<td>Date of:</td>
<td>2009-2020</td>
</tr>
</table>
</html>"));
  end Degradation;

  model Buffer
  "Source of substance bounded to constant amount of buffer to reach linear dependence between concentration and electrochemical potential"
    extends Icons.Buffer;
         extends Interfaces.PartialSubstanceInSolution(
                   a(start = a_start));

  parameter Modelica.Units.SI.MoleFraction a_start=1e-7
    "Initial value of mole fraction of the buffered substance";

  parameter Modelica.Units.SI.AmountOfSubstance BufferValue=0.001
    "Fixed buffer value (slope between amount of buffered substance and -log10(activity)) if useBufferValueInput=false"
    annotation (HideResult=true, Dialog(enable=not useBufferValueInput));

       parameter Boolean useBufferValueInput = false
    "Is buffer value of the substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        extends Interfaces.ConditionalKinetics(KC=1/(Modelica.Constants.R*298.15));

        Real bufferValue(final unit="1");

      Modelica.Blocks.Interfaces.RealInput bufferValueInput(
        final unit="mol/mol",
        start=BufferValue)=bufferValue
        if useBufferValueInput annotation (HideResult=true, Placement(transformation(
              extent={{-120,-20},{-80,20}})));

        Real xref;
  Modelica.Units.SI.AmountOfSubstance nFreeBuffer(start=-log10(a_start)
        *BufferValue) "amount of base molecules without H+";
  Modelica.Units.SI.MoleFraction xFreeBuffer;

  protected
  Modelica.Units.SI.MolarEnthalpy streamEnthalpy;

      constant Real InvLog_10=1/log(10);
  initial equation
      xFreeBuffer = -log10(a_start)*(bufferValue/solution.n);

  equation
      if not useBufferValueInput then
        bufferValue = BufferValue;
      end if;

      der(nFreeBuffer) = -port_a.q;
      // <- This is mathematically the same as two following lines. However, the differential solvers can handle the log10n much better. :-)
      //der(log10nFreeBuffer)=(InvLog_10)*(port_a.q/nFreeBuffer);
      //nFreeBuffer = 10^log10nFreeBuffer;

      xFreeBuffer = nFreeBuffer/solution.n;
     // port_a.q = (solution.n*KC)*(xFreeBuffer - xref);
      port_a.q = KC*(Modelica.Constants.R*solution.T*log(xFreeBuffer) - Modelica.Constants.R*solution.T*log(xref)); //alternative kinetics
      xref = -log10(a)*(bufferValue/solution.n);

    //solution flows
    streamEnthalpy = actualStream(port_a.h_outflow);

    solution.dH =streamEnthalpy*port_a.q - der(molarEnthalpy)*nFreeBuffer;
    solution.i = Modelica.Constants.F * z * port_a.q - Modelica.Constants.F*der(z)*nFreeBuffer;
    solution.dV = molarVolume * port_a.q - der(molarVolume)*nFreeBuffer;

    //extensive properties
    solution.nj=0;
    solution.mj=-nFreeBuffer*stateOfMatter.molarMassOfBaseMolecule(substanceData);
    solution.Vj=-nFreeBuffer*molarVolume;
    solution.Gj=-nFreeBuffer*port_a.u;
    solution.Qj=-Modelica.Constants.F*nFreeBuffer*z;
    solution.Ij=-(1/2) * ( nFreeBuffer * z^2);

      annotation ( Icon(coordinateSystem(
              preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
            graphics={
            Text(
              extent={{-82,62},{92,24}},
              textString="%name",
              lineColor={128,0,255})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Buffer;

  model SubstanceInflowT
    "Molar pump of substance at defined temperature to system"
    extends Interfaces.ConditionalSubstanceFlow;

  Interfaces.SubstancePort_b port_b "Outflow"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

   outer Modelica.Fluid.System system "System wide properties";

   replaceable package stateOfMatter =
      Chemical.Interfaces.Incompressible constrainedby Chemical.Interfaces.StateOfMatter
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

   parameter Modelica.Units.SI.Temperature T=system.T_ambient;
  equation
    port_b.q = -q;
    port_b.h_outflow = stateOfMatter.molarEnthalpy(substanceData,T);

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics={
          Rectangle(
            extent={{-100,-42},{100,40}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-48,20},{50,0},{-48,-21},{-48,20}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-220,-80},{220,-60}},
            textString="%name",
            lineColor={128,0,255})}),      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end SubstanceInflowT;
end Sources;
