within Chemical;
package Sensors "Sensors package for undirected chemical simulation"
  extends Modelica.Icons.SensorsPackage;

  import Chemical.Sensors.Internal;

  model SingleSensorSelect "Sensor with selectable measured quantity"
    extends Internal.PartialSensor;
    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    parameter Chemical.Sensors.Internal.Types.Quantities quantity "Quantity to be measured";
    parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimum density" annotation (Dialog(tab="Advanced", group="Regularization"));
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real value_0(unit=Chemical.Sensors.Internal.getUnit(quantity)) = 0 "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

    Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getUnit(quantity)) = value if outputValue "Measured quantity [variable]"
      annotation (Placement(
          transformation(extent={{80,-20},{120,20}}),
            iconTransformation(extent={{80,-20},{120,20}})));

    function getQuantity = Chemical.Sensors.Internal.getQuantity (
      redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
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

    direct_value = getQuantity(state, rear.r, quantity, rho_min);

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

  model TwoPhaseSensorSelect "Sensor for a selectable quantity of a twoPhaseMedium"
    extends Internal.PartialSensor(redeclare package Medium=Medium2Phase);

    import Quantities=Chemical.Sensors.Internal.Types.TwoPhaseQuantities;
    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    replaceable package Medium2Phase =
        Chemical.Media.myMedia.Interfaces.PartialTwoPhaseMedium
                                                       "Medium model"
      annotation (choicesAllMatching=true,
        Documentation(info="<html>
<u>Replaceable medium package for the sensor. Medium must be a TwoPase Medium.</u>
</html>"));

    parameter Quantities quantity "Quantity the sensor measures"
      annotation(choicesAllMatching=true);
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(choicesAllMatching=true, Dialog(tab="Initialization", enable=filter_output));
    parameter Real value_0(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity)) = 0 "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

    Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity)) = value if outputValue "Measured quantity [variable]"
      annotation (Placement(
          transformation(extent={{80,-20},{120,20}}),
            iconTransformation(extent={{80,-20},{120,20}})));

    Real value(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity));

  protected
    Real direct_value(unit=Chemical.Sensors.Internal.getTwoPhaseUnit(quantity));

    function getQuantity = Chemical.Sensors.Internal.getTwoPhaseQuantity (
      redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
      annotation (Documentation(info="<html>
    <u>This function computes the selected two-phase quantity from state.</u>
      </html>"));

  initial equation
    if filter_output and init==InitMode.steadyState then
      value= direct_value;
    elseif filter_output then
      value = value_0;
    end if;

  equation
    direct_value = getQuantity(state, quantity);

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
            points={{-100,-80},{100,-80}},
            color={158,66,200},
            thickness=0.5),
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
            extent={{0,19},{60,69}},
            textColor={175,175,175},
            textString="%quantity")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Undirected sensor for a vapor quality. It must be separate from SingleSensorSelect, because it needs a TwoPhaseMedium.</u>
</html>"));
  end TwoPhaseSensorSelect;

  model SensorState "Sensor for whole state"
    extends Internal.PartialSensor;

    replaceable package Medium =
        Chemical.Media.myMedia.Interfaces.PartialMedium
      "Medium model"
      annotation (choicesAllMatching=true,
        Documentation(info="<html>
        <u>Medium Model for the sensor. Make sure it is the same as for all lines the sensors input is connected.</u>
        </html>"));

    Chemical.Interfaces.StateOutput state_out(redeclare package stateOfMatter =
          stateOfMatter)                                                                       "Measured value [variable]"
      annotation (Placement(transformation(extent={{80,-20},{120,20}})));

  equation

    state_out.state = state;

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

  model SingleSensorX "Sensor for mass fraction of mixture"
    extends Internal.PartialSensor;

    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    replaceable package Medium =
        Chemical.Media.myMedia.Interfaces.PartialMedium
      "Medium model"
      annotation (choicesAllMatching=true,
        Documentation(info="<html>
        <u>Medium Model for the sensor. Make sure it is the same as for all lines the sensors input is connected.</u>
        </html>"));

    parameter Integer digits(min=0) = 1 "Number of displayed digits";
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real[Medium.nX] value_0(each unit="kg/kg") = Medium.X_default "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));
    parameter Integer row(min=1, max=Medium.nX) = 1 "Row of mass fraction vector to display";

    Modelica.Blocks.Interfaces.RealOutput value_out[Medium.nX](each unit="kg/kg") = value if outputValue "Measured value [variable]"
      annotation (Placement(transformation(extent={{80,-20},{120,20}})));

    output Real value[Medium.nX](each unit="kg/kg") "Computed value of the selected Quantity";
    output Real display_value(unit="kg/kg") = value[row] "Row of the value vector to display";

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Real direct_value[Medium.nX](each unit="kg/kg");

    function mfk = Chemical.Utilities.Functions.massFractionK (
                                                     redeclare package stateOfMatter =
            stateOfMatter);

  initial equation
    if filter_output and init==InitMode.steadyState then
      value= direct_value;
    elseif filter_output then
      value = value_0;
    end if;

  equation
    //OM fix
    //direct_value[1:Medium.nXi] = Medium.massFraction(state);

    for i in 1:Medium.nXi loop
      direct_value[i] = mfk(state, i);
    end for;

    if Medium.reducedX then
      direct_value[end] = 1-sum(direct_value[1:Medium.nXi]);
    end if;

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
          Rectangle(
            extent={{-60,30},{60,-30}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-60,30},{60,-30}},
            textColor={158,66,200},
            textString=DynamicSelect("value", String(display_value, format="1."+String(digits)+"f"))),
          Text(
            extent={{-26,22},{60,69}},
            textColor={175,175,175},
            textString="%row. mass-fraction"),
          Ellipse(
            extent={{-5,-75},{5,-85}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Sensor for measuring mass fraction X. Which row from X to display can be selected by the row parameter.</u>
<u>This sensor can be connected to a fluid stream without a junction.</u>
</html>"));
  end SingleSensorX;

  model SingleFlowSensor "Sensor for a selectable quantity associated with the massflow"
     extends Internal.PartialSensor;
    import Quantities=Chemical.Sensors.Internal.Types.MassFlowQuantities;
    import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;

    parameter Quantities quantity "Quantity the sensor measures";
    parameter Modelica.Units.SI.Density rho_min=dropOfCommons.rho_min "Minimum Density" annotation (Dialog(tab="Advanced", group="Regularization"));
    parameter Boolean outputValue = false "Enable sensor-value output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=outputValue));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real value_0(unit=Chemical.Sensors.Internal.getFlowUnit(quantity)) = 0 "Initial output state of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

    Modelica.Blocks.Interfaces.RealOutput value_out(unit=Chemical.Sensors.Internal.getFlowUnit(quantity)) = value if outputValue "Measured quantity [variable]"
      annotation (Placement(transformation(extent={{80,-20},{120,20}}),
          iconTransformation(extent={{80,-20},{120,20}})));

    output Real value(unit=Chemical.Sensors.Internal.getFlowUnit(quantity));

  protected
    Real direct_value(unit=Chemical.Sensors.Internal.getFlowUnit(quantity));

    function getQuantity = Chemical.Sensors.Internal.getFlowQuantity (
      redeclare package stateOfMatter = stateOfMatter) "Quantity compute function"
      annotation (Documentation(info="<html>
        <u>This function computes the selected quantity from state and massflow. rho_min is neddet for the computation of v. </u>
        </html>"));

  initial equation
    if filter_output and init==InitMode.steadyState then
      value= direct_value;
    elseif filter_output then
      value = value_0;
    end if;

  equation
    direct_value = getQuantity(state, rear.n_flow, quantity, rho_min);

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
            points={{0,-26},{0,-80}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{-6,-74},{6,-86}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5),
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
            extent={{0,25},{60,75}},
            textColor={175,175,175},
            textString="%quantity")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>A undirected sensor measuring a selectable flow quantity associated with the massflow. For some quatities several units are available.</u>
</html>"));
  end SingleFlowSensor;

  model MultiSensor_Tpm "Undirected Sensor for Temperature, potential and mass-flow"
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
    parameter Chemical.Sensors.Internal.Types.MassFlowUnit massFlowUnit = "(kg/s)" "Unit for potential measurement and output"
      annotation(choicesAllMatching = true, Evaluate = true);
    parameter Boolean outputTemperature = false "Enable temperature output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean outputChemicalPotential = false "Enable potential output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean outputMassFlowRate = false "Enable massFlow output"
      annotation(Dialog(group="Output Value"));
    parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
      annotation(Dialog(group="Output Value", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate)));
    parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
      annotation(Dialog(tab="Initialization", enable=filter_output));
    parameter Real u_0(final quantity="ChemicalPotential", final unit=potentialUnit) = 0 "Initial output potential of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Real T_0(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = 0 "Initial output Temperature of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Real n_flow_0(final quantity="MassFlowRate", final unit=massFlowUnit) = 0 "Initial output massflow of sensor"
      annotation(Dialog(tab="Initialization", enable=filter_output and init==InitMode.state));
    parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant"
      annotation (Dialog(tab="Advanced", enable=(outputTemperature or outputChemicalPotential or outputMassFlowRate) and filter_output));

    Modelica.Blocks.Interfaces.RealOutput T_out(final quantity="ThermodynamicTemperature", final unit=temperatureUnit) = T if outputTemperature "Measured temperature [variable]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,60}),
          iconTransformation(extent={{80,60},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput u_out(final quantity="ChemicalPotential", final unit=potentialUnit) = u if outputChemicalPotential "Measured potential [variable]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,20}),
          iconTransformation(extent={{80,0},{120,40}})));
    Modelica.Blocks.Interfaces.RealOutput n_flow_out(unit="kg/s") = n_flow if outputMassFlowRate
      "Measured mass-flow [kg/s]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-60}),
          iconTransformation(extent={{80,-60},{120,-20}})));

    output Real u(final quantity="ChemicalPotential", final unit=potentialUnit);
    output Real T(final quantity="ThermodynamicTemperature", final unit=temperatureUnit);
    output Real n_flow(final quantity="MassFlowRate", final unit=massFlowUnit);

  protected
    Real direct_p; // unit intentionally not set to avoid Warning
    Real direct_T; // unit intentionally not set to avoid Warning
    Real direct_n_flow; // unit intentionally not set to avoid Warning

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
      direct_T =Medium.temperature(state);
    elseif temperatureUnit == "degC" then
      direct_T =Modelica.Units.Conversions.to_degC(Medium.temperature(state));
    end if;

    if potentialUnit == "J/mol" then
      direct_p =state.u;
    elseif potentialUnit == "bar" then
      direct_p =Modelica.Units.Conversions.to_bar(state);
    end if;

    if massFlowUnit == "(kg/s)" then
        direct_n_flow = rear.n_flow;
    elseif massFlowUnit == "(g/s)" then
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
            textString=DynamicSelect("m", String(
                    n_flow,
                    format="1."+String(digits)+"f"))),
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
            textString="%massFlowUnit")}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>Undirected&nbsp;sensor&nbsp;for&nbsp;temperature,&nbsp;potential&nbsp;and&nbsp;mass-flow. Units can be selected.</u>
</html>"));
  end MultiSensor_Tpm;

  model UnidirectionalSensorAdapter "Adapter to connect a unidirectional sensor"
     extends Internal.PartialSensor;

    Chemical.Interfaces.Outlet outlet(redeclare package Medium=Medium)
      annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,40}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-40})));

  equation
    outlet.r = rear.r;
    outlet.state = state;

    annotation (Icon(
      graphics={
        Line(
          points={{0,-80},{0,-40}},
          color={158,66,200},
          thickness=0.5),
        Ellipse(
          extent={{-5,-85},{5,-75}},
          lineColor={158,66,200},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5)},
      coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,-20}})),
                                                    Diagram(
       coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,-20}})),
      Documentation(info="<html>
<u>A adapter to outputs the relevant state of the undirected flow, with r=0 at the outlet. It can be used to connect a unidirectional sensor to a undirected network.</u>
</html>"));
  end UnidirectionalSensorAdapter;

  package Tests "Tests package for the undirected sensor package"
  extends Modelica.Icons.ExamplesPackage;

    model TestSensors "Test for the undirected sensors"
      extends Modelica.Icons.Example;

      replaceable package Medium = Chemical.Media.myMedia.Water.StandardWater constrainedby
        Chemical.Media.myMedia.Interfaces.PartialTwoPhaseMedium
        "Medium model"
        annotation (Documentation(info="<html>
<u>Replaceable package with the medium model. Due to the vaporQuality sensor it must be a TwoPhaseMedium.</u>
</html>"));

      replaceable package Medium2 =
          Chemical.Media.myMedia.IdealGases.MixtureGases.CombustionAir                           constrainedby
        Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medium model"
        annotation (Documentation(info="<html>
<u>Replaceable package with the medium model. Due to the vaporQuality sensor it must be a TwoPhaseMedium.</u>
</html>"));

      Chemical.Processes.FlowResistance flowResistance(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=100,
        redeclare function pLoss = Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss (material=Chemical.Processes.Internal.Material.steel))
        annotation (Placement(transformation(extent={{-40,-8},{-20,12}})));
      Chemical.Boundaries.BoundaryRear boundary_rear(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=373.15,
        u0_par=100000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-120,2})));
      Chemical.Boundaries.BoundaryFore boundary_fore(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=303.15,
        u0_par=110000) annotation (Placement(transformation(extent={{84,-8},{104,12}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-130,22},{-110,42}})));
      Modelica.Blocks.Sources.Step step(
        height=-80000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{120,-4},{108,8}})));
      MultiSensor_Tpm multiSensor_Tpm(redeclare package stateOfMatter =
            stateOfMatter,
        temperatureUnit="degC",
        potentialUnit="bar",
        outputTemperature=false)
        annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
      SingleSensorSelect singleSensorSelect(redeclare package stateOfMatter =
            stateOfMatter,
          quantity=Chemical.Sensors.Internal.Types.Quantities.u_bar)
        annotation (Placement(transformation(extent={{-10,0},{10,20}})));
      UnidirectionalSensorAdapter unidirectionalSensorAdapter(
        redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(extent={{20,0},{40,8}})));
      Chemical.Sensors.TwoPhaseSensorSelect sensor_vaporQuality1(
        redeclare package stateOfMatter = stateOfMatter, quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg)
        annotation (Placement(transformation(extent={{50,12},{70,32}})));
      SingleFlowSensor singleFlowSensor(redeclare package stateOfMatter =
            stateOfMatter,                                                               quantity=Chemical.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps)
        annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      TwoPhaseSensorSelect sensor_vaporQuality2(
        redeclare package Medium2Phase = Medium,
        quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg,
        redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(extent={{50,0},{70,20}})));
      Chemical.Processes.FlowResistance flowResistance1(
        redeclare package stateOfMatter = stateOfMatter,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=100,
        redeclare function pLoss = Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss (material=Chemical.Processes.Internal.Material.steel))
        annotation (Placement(transformation(extent={{-40,52},{-20,72}})));
      Chemical.Boundaries.BoundaryRear boundary_rear1(
        redeclare package stateOfMatter = stateOfMatter,
        T0_par=373.15,
        u0_par=100000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-120,62})));
      Chemical.Boundaries.BoundaryFore boundary_fore1(
        redeclare package stateOfMatter = stateOfMatter,
        potentialFromInput=true,
        T0_par=303.15,
        u0_par=110000) annotation (Placement(transformation(extent={{84,52},{104,72}})));
      Modelica.Blocks.Sources.Step step1(
        height=-80000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{120,56},{108,68}})));
      MultiSensor_Tpm multiSensor_Tpm1(
        redeclare package stateOfMatter = stateOfMatter,
        temperatureUnit="degC",
        potentialUnit="bar",
        outputTemperature=true,
        filter_output=true)
        annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
      SingleSensorSelect singleSensorSelect1(
        redeclare package stateOfMatter = stateOfMatter,
        outputValue=true,
        quantity=Chemical.Sensors.Internal.Types.Quantities.u_bar,
        filter_output=true,
        init=Chemical.Sensors.Internal.Types.InitializationModelSensor.state,
        value_0=1) annotation (Placement(transformation(extent={{-10,60},{10,80}})));
      SingleFlowSensor singleFlowSensor1(
        redeclare package stateOfMatter = stateOfMatter,
        outputValue=true,
        quantity=Chemical.Sensors.Internal.Types.MassFlowQuantities.H_flow_Jps,
        filter_output=true) annotation (Placement(transformation(extent={{-70,60},{-50,80}})));
      TwoPhaseSensorSelect sensor_vaporQuality4(
        redeclare package Medium2Phase = Medium,
        outputValue=true,
        quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg,
        filter_output=true,
        redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(extent={{50,60},{70,80}})));
      UnidirectionalSensorAdapter unidirectionalSensorAdapter1(
        redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(extent={{20,64},{40,56}})));
      Chemical.Sensors.DifferenceTwoPhaseSensorSensorSelect differenceTwoPhaseSensorSensorSelect(
        redeclare package MediumA = Medium,
        redeclare package MediumB = Medium,
        quantity=Chemical.Sensors.Internal.Types.TwoPhaseQuantities.u_sat_Pa,
        outputValue=true,
        filter_output=true) annotation (Placement(transformation(extent={{50,52},{70,32}})));
      Chemical.Boundaries.BoundaryRear boundary_rear2(
        redeclare package Medium = Medium2,
        T0_par=373.15,
        u0_par=200000,
        Xi0_par={0.2,0.8}) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-120,-30})));
      Chemical.Boundaries.BoundaryFore boundary_fore2(
        redeclare package Medium = Medium2,
        potentialFromInput=false,
        T0_par=303.15,
        u0_par=100000) annotation (Placement(transformation(extent={{86,-40},{106,-20}})));
      Chemical.Processes.FlowResistance flowResistance2(
        redeclare package Medium = Medium2,
        initM_flow=Chemical.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=100,
        redeclare function pLoss = Chemical.Processes.Internal.FlowResistance.laminarTurbulentChemicalPotentialLoss (material=Chemical.Processes.Internal.Material.steel))
        annotation (Placement(transformation(extent={{-38,-40},{-18,-20}})));
      SingleSensorX singleSensorX(redeclare package Medium = Medium2) annotation (Placement(transformation(extent={{-100,-32},{-80,-12}})));
      SingleSensorX singleSensorX1(
        redeclare package Medium = Medium2,
        digits=2,
        outputValue=true,
        filter_output=true,
        row=2) annotation (Placement(transformation(extent={{-70,-32},{-50,-12}})));
      SensorState sensorState(redeclare package Medium = Medium2) annotation (Placement(transformation(extent={{20,-32},{40,-12}})));
    equation
      connect(step.y, boundary_fore.u0_var)
        annotation (Line(points={{107.4,2},{102,2},{102,8},{96,8}},
                                                       color={0,0,127}));
      connect(singleSensorSelect.rear, flowResistance.fore) annotation (Line(
          points={{-10,2},{-20,2}},
          color={158,66,200},
          thickness=0.5));
      connect(singleSensorSelect.fore, unidirectionalSensorAdapter.rear)
        annotation (Line(
          points={{10,2},{20,2}},
          color={158,66,200},
          thickness=0.5));
      connect(sensor_vaporQuality1.inlet, unidirectionalSensorAdapter.outlet)
        annotation (Line(
          points={{50,22},{30,22},{30,6}},
          color={158,66,200},
          thickness=0.5));
      connect(multiSensor_Tpm.fore, singleFlowSensor.rear) annotation (Line(
          points={{-80,2},{-70,2}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance.rear, singleFlowSensor.fore) annotation (Line(
          points={{-40,2},{-50,2}},
          color={158,66,200},
          thickness=0.5));
      connect(unidirectionalSensorAdapter.fore, sensor_vaporQuality2.rear)
        annotation (Line(
          points={{40,2},{50,2}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore.rear, sensor_vaporQuality2.fore)
        annotation (Line(
          points={{84,2},{70,2}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear.fore, multiSensor_Tpm.rear) annotation (Line(
          points={{-110,2},{-100,2}},
          color={158,66,200},
          thickness=0.5));
      connect(step1.y, boundary_fore1.u0_var) annotation (Line(points={{107.4,62},{102,62},{102,68},{96,68}},
                                                                                              color={0,0,127}));
      connect(singleSensorSelect1.rear, flowResistance1.fore)
        annotation (Line(
          points={{-10,62},{-20,62}},
          color={158,66,200},
          thickness=0.5));
      connect(multiSensor_Tpm1.fore, singleFlowSensor1.rear)
        annotation (Line(
          points={{-80,62},{-70,62}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance1.rear, singleFlowSensor1.fore)
        annotation (Line(
          points={{-40,62},{-50,62}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_fore1.rear, sensor_vaporQuality4.fore)
        annotation (Line(
          points={{84,62},{70,62}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear1.fore, multiSensor_Tpm1.rear)
        annotation (Line(
          points={{-110,62},{-100,62}},
          color={158,66,200},
          thickness=0.5));
      connect(singleSensorSelect1.fore, unidirectionalSensorAdapter1.rear)
        annotation (Line(
          points={{10,62},{20,62}},
          color={158,66,200},
          thickness=0.5));
      connect(unidirectionalSensorAdapter1.fore, sensor_vaporQuality4.rear)
        annotation (Line(
          points={{40,62},{50,62}},
          color={158,66,200},
          thickness=0.5));
      connect(differenceTwoPhaseSensorSensorSelect.inletA, unidirectionalSensorAdapter.outlet)
        annotation (Line(
          points={{50.4,38},{30,38},{30,6}},
          color={158,66,200},
          thickness=0.5));
      connect(differenceTwoPhaseSensorSensorSelect.inletB, unidirectionalSensorAdapter1.outlet)
        annotation (Line(
          points={{50.4,46},{30,46},{30,58}},
          color={158,66,200},
          thickness=0.5));
      connect(singleSensorX.rear, boundary_rear2.fore) annotation (Line(
          points={{-100,-30},{-110,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance2.rear, singleSensorX1.fore) annotation (Line(
          points={{-38,-30},{-50,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(singleSensorX1.rear, singleSensorX.fore) annotation (Line(
          points={{-70,-30},{-80,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(flowResistance2.fore, sensorState.rear) annotation (Line(
          points={{-18,-30},{20,-30}},
          color={158,66,200},
          thickness=0.5));
      connect(sensorState.fore, boundary_fore2.rear) annotation (Line(
          points={{40,-30},{86,-30}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-60},{140,80}})),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Test&nbsp;for&nbsp;the&nbsp;undirected&nbsp;sensors.</u>
<u>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></u>
</html>"));
    end TestSensors;
    annotation (Documentation(info="<html>
<u>Tests package for the undirected sensor package.</u>
</html>"));
  end Tests;

  package Internal "Partials and functions"
    extends Modelica.Icons.InternalPackage;

    partial model PartialSensor "Partial undirected sensor"
      replaceable package Medium =
          Chemical.Media.myMedia.Interfaces.PartialMedium
        "Medium model" annotation (choicesAllMatching=true, Documentation(
            info="<html>
<u>Replaceable medium package for the sensor.</u>
</html>"));

      parameter Integer digits(min=0) = 1 "Number of displayed digits"
        annotation(Dialog(group="Sensor display"));
      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation (Dialog(tab="Advanced", group="Regularization"));

      Chemical.Interfaces.Rear rear(redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={-100,-80})));
      Chemical.Interfaces.Fore fore(redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(extent={{-20,-20},{20,20}}, origin={100,-80})));

    /*  function regStepSt = Undirected.Internal.regStepState (
    redeclare package stateOfMatter = stateOfMatter) "RegStep function for a state"
    annotation (Documentation(info="<html>
<u><span style=\"font-family: Courier New;\">RegStep function for a state. The medium of the sensor is used and given to the function.</span></u>
</html>"));
*/
      Medium.ThermodynamicState state = Chemical.Interfaces.SubstanceState(u=u_reg,h= h_reg); //= regStepSt(rear.n_flow, rear.state_forwards, rear.state_rearwards, n_flow_reg);

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
      for i in 1:Medium.nXi loop
        Xi_reg[i] =Chemical.Utilities.Internal.regStep(
              rear.n_flow,
              Xi_forwards[i],
              Xi_rearwards[i],
              n_flow_reg);
      end for;

      fore.state_forwards = rear.state_forwards;
      rear.state_rearwards = fore.state_rearwards;
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
