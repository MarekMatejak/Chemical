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

  model MoleFractionSensor "Measure of mole fraction"
    extends Modelica.Icons.RoundSensor;
    extends Internal.SubstanceSensor;

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

  model MolalitySensor "Measure of molality of the substance"
    extends Modelica.Icons.RoundSensor;
    extends Internal.SubstanceSensor;

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

  equation

    molality = x * (inlet.solution.n/inlet.solution.m);

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
    extends Internal.SubstanceSensor;


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

  equation

    molarConcentration = x * (inlet.solution.n/inlet.solution.V);

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
    extends Internal.SubstanceSensor;


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

    massFraction = (x / stateOfMatter.specificAmountOfParticles(substanceData)) * (inlet.solution.n/inlet.solution.m);

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
    extends Internal.SubstanceSensor;

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

    partialPressure = x*inlet.solution.p;

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

  package Internal

    model SubstanceSensor "Base class for sensor based on inlet substance and solution properties"
      extends Interfaces.PartialSubstanceInlet;
    equation
      inlet.n_flow=0;
    end SubstanceSensor;

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
      input Modelica.Units.SI.ChemicalPotential r "Inertial electro-chemical potential";
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
      Modelica.Units.SI.ChemicalPotential u_Pure;
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


      u_Pure := stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength)
       + z*Modelica.Constants.F*electricPotential;

      a := exp((u - u_Pure)/(Modelica.Constants.R*temperature));
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
    extends Internal.SubstanceSensor;

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

    function getQuantity = Internal.getQuantity (redeclare package stateOfMatter =
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

    direct_value = getQuantity(inlet.state.u, inlet.state.h, inlet.r, quantity, substanceData, temperature, pressure, electricPotential, moleFractionBasedIonicStrength,
     inlet.solution.m, inlet.solution.n, inlet.solution.V);


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
