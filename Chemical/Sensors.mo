within Chemical;
package Sensors "Sensors package for undirected chemical simulation"
  extends Modelica.Icons.SensorsPackage;

  import Chemical.Sensors.Internal;

  model SingleSensorSelect "Sensor with selectable measured quantity"
    extends Internal.PartialSensor;
    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    parameter Chemical.Sensors.Internal.Types.Quantities quantity "Quantity to be measured";
  //  parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimum density" annotation (Dialog(tab="Advanced", group="Regularization"));
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real value_0(unit=Chemical.Sensors.Internal.getUnit(quantity)) = 0 "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=dropOfCommons.TC "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

    Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getUnit(quantity)) = value if outputValue "Measured quantity [variable]"
      annotation (Placement(
          transformation(extent={{80,-20},{120,20}}),
            iconTransformation(extent={{80,-20},{120,20}})));

    function getQuantity = Chemical.Sensors.Internal.getQuantity "Quantity compute function"
      annotation (Documentation(info="<html>
      <u>This function computes the selected quantity from state. r and rho_min are neddet for the quantities r/u_total and v respectively.</u>
      </html>"));

    Real value(unit=Chemical.Sensors.Internal.getUnit(quantity));

  protected
    Real direct_value(unit=Chemical.Sensors.Internal.getUnit(quantity));

  initial equation
    if filter_output and init==InitMode.steadyState then
      value= direct_value;
    elseif filter_output then
      value = value_0;
    end if;

  equation

    /*
    input Chemical.Interfaces.SubstanceState substanceState "Substance state";
  input Modelica.Units.SI.ChemicalPotential r "Inertial electro-chemical potential";
  input Quantities quantity "What to measure?";
  input Chemical.Interfaces.Definition substance "Substance definition";
  input Chemical.Interfaces.SolutionState solution "Schemica solution state";
  */
    direct_value = getQuantity(state, rear.r, quantity, rear.definition, rear.solution_forwards);

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
            points={{0,0},{0,-80}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-5,-75},{5,-85}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-60,30},{60,-30}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-60,30},{60,-30}},
            textColor={158,66,200},
            textString=DynamicSelect("value", String(value, format="1."+String(digits)+"f"))),
          Text(
            extent={{2,17},{62,67}},
            textColor={175,175,175},
            textString="%quantity")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Undirected sensor for a single selectable quantity. For some quatities several units are available.</u>
</html>"));
  end SingleSensorSelect;

  model SensorState "Sensor for whole state"
    extends Internal.PartialSensor;


    Chemical.Interfaces.SubstanceStateOutput state_out   "Measured value [variable]"
      annotation (Placement(transformation(extent={{80,-20},{120,20}})));

  equation

    state_out = state;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-54,24},{66,-36}},
            lineColor={0,0,0},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{0,0},{0,-80}},
            color={158,66,200},
            thickness=0.5),
          Rectangle(
            extent={{-60,30},{60,-30}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-60,30},{60,-30}},
            textColor={158,66,200},
            textString="state"),
          Ellipse(
            extent={{-5,-75},{5,-85}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Sensor for measuring the full state.</u>
<u>This sensor can be connected to a fluid stream without a junction.</u>
</html>"));
  end SensorState;

  model MultiSensor_Tpn "Undirected Sensor for Temperature, potential and molar-flow"
    extends Internal.PartialSensor;

    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    parameter Chemical.Sensors.Internal.Types.TemperatureUnit temperatureUnit="K" "Unit for the temperature output"
      annotation (
      Dialog(group="Units"),
      choicesAllMatching=true,
      Evaluate=true);
    parameter Chemical.Sensors.Internal.Types.ChemicalPotentialUnit potentialUnit="J/mol" "Unit for the potential output"
      annotation (
      Dialog(group="Units"),
      choicesAllMatching=true,
      Evaluate=true);
    parameter Chemical.Sensors.Internal.Types.FlowUnit flowUnit = "(mol/s)" "Unit for potential measurement and output"
      annotation(choicesAllMatching = true, Evaluate = true);
    parameter Boolean outputTemperature = false "Enable temperature output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean outputChemicalPotential = false "Enable potential output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean outputMolarFlowRate = false "Enable molarFlow output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate)));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real u_0(final quantity="ChemicalPotential", final unit=potentialUnit) = 0 "Initial output potential of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Real T_0(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = 0 "Initial output Temperature of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Real n_flow_0(final quantity="MolarFlowRate", final unit=flowUnit) = 0 "Initial output molarflow of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=dropOfCommons.TC "PT1 time constant"
      annotation (Dialog(tab="Advanced", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate) and filter_output));

    Modelica.Blocks.Interfaces.RealOutput T_out(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = T if outputTemperature "Measured temperature [variable]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,60}),
          iconTransformation(extent={{80,60},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput u_out(final quantity="ChemicalPotential", final unit=potentialUnit) = u if outputChemicalPotential "Measured potential [variable]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,20}),
          iconTransformation(extent={{80,0},{120,40}})));
    Modelica.Blocks.Interfaces.RealOutput n_flow_out(unit="kg/s") = n_flow if outputMolarFlowRate
      "Measured molar-flow [molar/s]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-60}),
          iconTransformation(extent={{80,-60},{120,-20}})));

    output Real u(final quantity="ChemicalPotential", final unit=potentialUnit);
    output Real T(final quantity="ThermodynamicTemperature", final unit=temperatureUnit);
    output Real n_flow(final quantity="MassFlowRate", final unit=flowUnit);

  protected
    Real direct_p; // unit intentionally not set to avoid Warning
    Real direct_T; // unit intentionally not set to avoid Warning
    Real direct_n_flow; // unit intentionally not set to avoid Warning

    Chemical.Interfaces.SolutionState solution = rear.solution_forwards;

  initial equation
    if filter_output and init==InitMode.steadyState then
      u=direct_p;
      T=direct_T;
      n_flow=direct_n_flow;
    elseif filter_output then
      u=u_0;
      T=T_0;
      n_flow=n_flow_0;
    end if;

  equation
    if temperatureUnit == "K" then
      direct_T = solution.T;
    elseif temperatureUnit == "degC" then
      direct_T =Modelica.Units.Conversions.to_degC(solution.T);
    end if;

    if potentialUnit == "J/mol" then
      direct_p =state.u;
    elseif potentialUnit == "bar" then
      direct_p =Modelica.Units.Conversions.to_bar(solution.p);
    end if;

    if flowUnit == "(mol/s)" then
        direct_n_flow = rear.n_flow;
    elseif flowUnit == "(mmol/s)" then
        direct_n_flow = rear.n_flow*1000;
    end if;

    if filter_output then
      der(u) * TC = direct_p-u;
      der(T) * TC = direct_T-T;
      der(n_flow) * TC = direct_n_flow-n_flow;
    else
      u = direct_p;
      T = direct_T;
      n_flow = direct_n_flow;
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-54,94},{66,-66}},
            lineColor={0,0,0},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-60,100},{60,-60}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(points={{0,-60},{0,-80}}, color={0,0,0}),
          Ellipse(
            extent={{-6,-74},{6,-86}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
          Text(
            extent={{-60,100},{60,50}},
            textColor={158,66,200},
            textString=DynamicSelect("T", String(
                    T,
                    format="1."+String(digits)+"f"))),
          Text(
            extent={{-60,50},{60,0}},
            textColor={158,66,200},
            textString=DynamicSelect("u", String(
                    u,
                    format="1."+String(digits)+"f"))),
          Text(
            extent={{-60,0},{60,-50}},
            textColor={158,66,200},
            textString="n"),
          Text(
            extent={{-120,100},{-60,48}},
            textColor={175,175,175},
            textString="%temperatureUnit"),
          Text(
            extent={{-120,52},{-60,0}},
            textColor={175,175,175},
            textString="%potentialUnit"),
          Text(
            extent={{-120,0},{-60,-52}},
            textColor={175,175,175},
            textString="%flowUnit")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Undirected&nbsp;sensor&nbsp;for&nbsp;temperature,&nbsp;potential&nbsp;and&nbsp;mass-flow. Units can be selected.</u>
</html>"));
  end MultiSensor_Tpn;

  package Internal "Partials and functions"
    extends Modelica.Icons.InternalPackage;

    partial model PartialSensor "Partial undirected sensor"


      parameter Integer digits(min=0) = 1 "Number of displayed digits"
        annotation(Dialog(group="Sensor display"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation (Dialog(tab="Advanced", group="Regularization"));

      Chemical.Interfaces.Rear rear
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,-80})));
      Chemical.Interfaces.Fore fore
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-80})));

    /*  function regStepSt = Undirected.Internal.regStepState (
    redeclare package stateOfMatter = stateOfMatter) "RegStep function for a state"
    annotation (Documentation(info="<html>
<u><span style=\"font-family: Courier New;\">RegStep function for a state. The medium of the sensor is used and given to the function.</span></u>
</html>"));
*/
      Chemical.Interfaces.SubstanceState state = Chemical.Interfaces.SubstanceState(u=u_reg,h= h_reg); //= regStepSt(rear.n_flow, rear.state_forwards, rear.state_rearwards, n_flow_reg);

    protected
      outer Chemical.DropOfCommons dropOfCommons;
     Modelica.Units.SI.ChemicalPotential u_reg=Chemical.Utilities.Internal.regStep(
                rear.n_flow,
                rear.state_forwards.u,
                rear.state_rearwards.u,
                n_flow_reg);
      Modelica.Units.SI.MolarEnthalpy h_reg=Chemical.Utilities.Internal.regStep(
                rear.n_flow,
                rear.state_forwards.h,
                rear.state_rearwards.h,
                n_flow_reg);

    equation

      fore.state_forwards = rear.state_forwards;
      rear.state_rearwards = fore.state_rearwards;
      connect(rear.definition,fore.definition);
      connect(rear.solution_forwards,fore.solution_forwards);
      connect(rear.solution_rearwards,fore.solution_rearwards);
      fore.r = rear.r;
      fore.n_flow + rear.n_flow = 0;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,-80},{100,-80}},
              color={158,66,200},
              thickness=0.5)}),
           Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<u>Parent class of all undirected sensors.</u>
</html>"));
    end PartialSensor;

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
      type TemperatureUnit
      extends String;

        annotation (choices(
            choice="degC" "°C",
            choice="K" "Kelvin"));

      end TemperatureUnit;

      type ChemicalPotentialUnit
      extends String;

        annotation (choices(
            choice="J/mol" "J/mol",
            choice="kJ/mol" "kJ/mol"));

      end ChemicalPotentialUnit;

      type FlowQuantities = enumeration(
          n_flow_molps "Molar flow (mol/s)",
          n_flow_mmolps "Molar flow (mmol/s)",
          m_flow_kgps "Mass flow (kg/s)",
          m_flow_gps "Mass flow (g/s)",
          H_flow_Jps "Enthalpy flow (J/s)",
          S_flow_JpKs "Enthopy flow (J/(K.s))",
          Cp_flow_JpKs "Heat capacity flow (J/(K.s))",
          V_flow_m3ps "Volume flow (m3/s)",
          V_flow_lpMin "Volume flow (l/min)");
      type FlowUnit
      extends String;

        annotation (choices(
            choice="(mol/s)" "mol/s",
            choice="(mmol/s)" "mmol/s"));

      end FlowUnit;
    end Types;

    function getQuantity "Computes selected quantity from state"
      extends Modelica.Icons.Function;

      import Chemical.Sensors.Internal.Types.Quantities;

      input Chemical.Interfaces.SubstanceState substanceState "Substance state";
      input Modelica.Units.SI.ChemicalPotential r "Inertial electro-chemical potential";
      input Quantities quantity "What to measure?";
      input Chemical.Interfaces.Definition substance "Substance definition";
      input Chemical.Interfaces.SolutionState solution "Schemica solution state";


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
      z :=  Chemical.Interfaces.Properties.chargeNumberOfIon( substance, solution);
      gamma := Chemical.Interfaces.Properties.activityCoefficient( substance, solution);

      u_Pure := Chemical.Interfaces.Properties.chemicalPotentialPure( substance, solution)
       + z*Modelica.Constants.F*solution.v;

      a := exp((substanceState.u - u_Pure)/(Modelica.Constants.R*solution.T));
      x := a/gamma;

      if quantity == Quantities.c_molpm3 then
        value := (x * solution.n)/solution.V;
      elseif quantity == Quantities.X_kgpkg then
        value := ((x * solution.n)/solution.m)/Chemical.Interfaces.Properties.specificAmountOfParticles( substance, solution);
      elseif quantity == Quantities.b_molpkg then
        value := (x * solution.n)/solution.m;
      elseif quantity == Quantities.x_molpmol then
        value := x;
      elseif quantity == Quantities.p_Pa then
        value := x*solution.p;
      elseif quantity == Quantities.p_mmHg then
        value := x*solution.p * (760/101325);
      elseif quantity == Quantities.p_bar then
        value := Modelica.Units.Conversions.to_bar(x*solution.p);
      elseif quantity == Quantities.u_Jpmol then
        value := substanceState.u;
      elseif quantity == Quantities.u_kJpmol then
        value := substanceState.u/1000;
      elseif quantity == Quantities.h_Jpmol then
        value := substanceState.h;
      elseif quantity == Quantities.s_JpmolK then
        value := (substanceState.h-substanceState.u)/solution.T;
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
    annotation (Documentation(info="<html>
<u>Partials and functions needet for undirected sensors.</u>
</html>"));
  end Internal;
  annotation (Documentation(revisions="<html>
<u><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</u>

</html>", info="<html>
<u>
Sensors package for undirected chemical simulation. For undirected flow, the massflow
must always be measured, since it determines if the forward or the backward flow is valid.
Therefore the fluid must flow through the sensor.
</u>
</html>"));
end Sensors;
