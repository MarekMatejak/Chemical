within Chemical;
package Boundaries "Boundary models for undirected chemical simulation"
  extends Modelica.Icons.SourcesPackage;

  model Substance "Substance"
    extends Icons.Substance;
    extends Internal.PartialSubstance(
      useFore=false,
      useSolution=false,
      useRear=false,
      final m_start=if preferMass then mass_start else
       amountOfSubstance_start*Chemical.Interfaces.Properties.molarMassOfBaseMolecule(substanceDefinition));

    import Chemical.Utilities.Types.InitializationSubstance;

    parameter Boolean preferMass=true  "prefere state as mass, otherwise amountOfSubstance"
      annotation (HideResult=true, Evaluate=true, choices(checkBox=true), Dialog(group="Substance"));

    parameter Modelica.Units.SI.Mass mass_start=dropOfCommons.DefaultMass "Initial mass of the substance"
      annotation (HideResult=true, Dialog(group="Substance", enable=preferMass));

    parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=dropOfCommons.DefaultAmount
    "Initial amount of substance base molecules"
      annotation (HideResult=true, Dialog(group="Substance", enable=(not preferMass)));

    Modelica.Units.SI.Mass mass=n*molarMassOfBaseMolecule "Mass";

    parameter InitializationSubstance initAmount=Chemical.Utilities.Types.InitializationSubstance.state "Initialization"
      annotation (HideResult=true, Dialog(tab="Initialization"));

    parameter Modelica.Units.SI.MolarFlowRate change_start=0
    "Initial molar change of substance base molecules"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(initAmount ==InitializationSubstance.derivative)));

  protected

    Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = Chemical.Interfaces.Properties.molarMassOfBaseMolecule(substanceDefinitionVar);

    Real logn(stateSelect=StateSelect.prefer, start=log(amountOfSubstance_start))
    "Natural logarithm of the amount of base molecules in solution";

    Real logm(stateSelect=StateSelect.prefer, start=log(mass_start))
    "Natural logarithm of the substance mass in solution";

  initial equation

    if initAmount ==InitializationSubstance.steadyStateForwards  then
      substance.u = state_in_rear.u;
    elseif initAmount ==InitializationSubstance.steadyStateRearwards  then
      substance.u = state_in_fore.u;
    elseif initAmount ==InitializationSubstance.state  and not preferMass then
      logn=log(amountOfSubstance_start);
    elseif initAmount ==InitializationSubstance.state  and preferMass then
      logm=log(mass_start);
    elseif initAmount ==InitializationSubstance.derivative  then
      n_flow = change_start;
    end if;

  equation

    //The main accumulation equation is "der(n)=n_flow"
    // However, the numerical solvers can handle it during equilibration of chemical potential in form of log(n) much better. :-)
    if preferMass then
      der(logm) = (n_flow/n) "accumulation of m=exp(logm) [kg]";
    else
      der(logn) = (n_flow/n) "accumulation of n=exp(logn) [mol]";
    end if;
    //der(n) = n_flow;
    n = exp(logn);
    n = exp(logm)*1/molarMassOfBaseMolecule;

       annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={128,0,255},
            textString="%name")}), Documentation(revisions="<html>
<p>2009-2025 by Marek Matejak, Ph.D. </p>
</html>",   info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solutionState. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
<p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
<p>x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
<p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
<p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of x will be the same as molality.</p>
<p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of x will be the same as molarity.</p>
<p><br><br>Definition of electro-chemical potential:</p>
<h4>u = u&deg; + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg; = DfG = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. relative eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg; .. chemical potential of pure substance</p>
<p>DfG .. free Gibbs energy of formation of the substance</p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
  end Substance;

  model ElectronTransfer
    extends Icons.ElectronTransfer;
    extends Chemical.Boundaries.Internal.PartialBoundaryBase;
    outer Modelica.Fluid.System system "System wide properties";

    import Chemical.Utilities.Types.InitializationSubstance;

    parameter InitializationSubstance init=Chemical.Utilities.Types.InitializationSubstance.state "Initialization"
      annotation (HideResult=true, Dialog(tab="Initialization"));

    parameter Modelica.Units.SI.ElectricCurrent i_start=0 "Initial electric current"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(init ==InitializationSubstance.derivative)));

    parameter Modelica.Units.SI.ElectricPotential v_start=0 "Initial electric potential"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(init ==InitializationSubstance.state)));

    parameter Boolean useRear = false "Use rearwards conector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));
    parameter Boolean useFore = true "Use forwards connector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));


    Chemical.Interfaces.Rear rear(
      n_flow=n_flow_rear,
      r=r_rear_port,
      state_rearwards=state_out,
      solution_rearwards=solutionState) if useRear annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

    Chemical.Interfaces.Fore fore(
      n_flow=n_flow_fore,
      r=r_fore_port,
      state_forwards=state_out,
      solution_forwards=solutionState,
      definition=definition) if useFore annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

     Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
        Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
            extent={{-10,88},{10,108}})));

     Chemical.Interfaces.SolutionPort solution(
        T=solutionState.T,
        p=solutionState.p,
        v=solutionState.v,
        n=solutionState.n,
        m=solutionState.m,
        V=solutionState.V,
        Q=solutionState.Q,
        I=solutionState.I,
        G=solutionState.G,
        dH=0,
        dV=0,
        nj=0,
        mj=0,
        Vj=0,
        Gj=0,
        Qj=0,
        Ij=0)
      "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

  protected

    parameter Chemical.Interfaces.Definition definition=Chemical.Substances.Solid.e "Definition of the substance"
      annotation (choicesAllMatching=true, Dialog(enable=not useRear));


    Chemical.Interfaces.SolutionState solutionState; // = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible);


  initial equation
     if init == InitializationSubstance.steadyStateForwards  then
      state_out.u = state_in_rear.u;
    elseif init ==InitializationSubstance.steadyStateRearwards  then
      state_out.u = state_in_fore.u;
    elseif init ==InitializationSubstance.state then
      state_out.u = definition.data.z*Modelica.Constants.F*v_start;
    elseif init ==InitializationSubstance.derivative  then
      n_flow = i_start/(definition.data.z*Modelica.Constants.F);
    end if;

  equation
     //electric
    pin.v = solution.v;
    pin.i + definition.data.z*Modelica.Constants.F*n_flow + solution.i = 0;

    state_out.u = definition.data.z*Modelica.Constants.F*solutionState.v;
    state_out.h = definition.data.z*Modelica.Constants.F*solutionState.v;


    connect(state_in_rear,rear.state_forwards);
    connect(state_in_fore,fore.state_rearwards);

    if not useRear then
      r_rear_port = 0;
      n_flow_rear = 0;
      state_in_rear.h = 0;
    end if;

    if not useFore then
      r_fore_port = 0;
      n_flow_fore = 0;
      state_in_fore.h = 0;
    end if;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ElectronTransfer;

  model ExternalSubstance "Constant source of molar concentration"
    extends Internal.PartialSubstance(
      final m_start=1,
      substance(final SolutionObserverOnly=true));

   /* parameter stateOfMatter.SubstanceDataParameters substanceData
 "Definition of the substance"
    annotation (choicesAllMatching = true, Dialog(enable=not useRear));
*/
    parameter Chemical.Boundaries.Internal.Types.ConcentrationQuantities quantity "Concentration quantity";

    parameter Real FixedValue
    "Fixed value of concentration in selected quantity if useVariableInput=false"
      annotation (HideResult=true, Dialog(enable=not useVariableInput));

    parameter Boolean useVariableInput = false
    "Is amount of substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput VariableInput=val if useVariableInput
      annotation (HideResult=true, Placement(transformation(extent={{-130,56},{-90,96}})));

    Real  value(unit=Chemical.Boundaries.Internal.getUnit(           quantity));
  equation

    if not useVariableInput then
      value=FixedValue;
    end if;

    if quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3 then
      value =  1*(substance.x * solutionState.n)/solutionState.V;
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.X_kgpkg then
      value = 1*((substance.x * solutionState.n)/solutionState.m)/Chemical.Interfaces.Properties.specificAmountOfParticles(substanceDefinitionVar,
     solutionState);
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.b_molpkg then
      value = 1* (substance.x * solutionState.n)/solutionState.m;
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.x_molpmol then
      value = substance.x;
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_Pa then
      value =  1*substance.x*solutionState.p;
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_kPa then
      value*1000 =  substance.x*solutionState.p;
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg then
      value =  substance.x*solutionState.p * (760/101325);
    elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_bar then
      value = 1* Modelica.Units.Conversions.to_bar(substance.x*solutionState.p);
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
  end ExternalSubstance;

  model ExternalGas "Gas substance with defined partial pressure"
    extends Internal.PartialSubstance(
      solutionParam=Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Gas),
      m_start=1,
      substance(final SolutionObserverOnly=true));


    parameter Boolean usePartialPressureInput = false
    "=true, if fixed partial pressure is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Pressure PartialPressure
      "Fixed partial pressure if usePartialPressureInput=false" annotation (
       HideResult=true, Dialog(enable=not usePartialPressureInput));

    Modelica.Blocks.Interfaces.RealInput partialPressure(start=
          PartialPressure, final unit="Pa")=p if usePartialPressureInput
    "Partial pressure of gas = total pressure * gas fraction"
      annotation (HideResult=true,Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-100,72}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-100,72})));

    Modelica.Units.SI.Pressure p "Current partial pressure";

  equation

    if not usePartialPressureInput then
      p=PartialPressure;
    end if;

    //mole fraction
    substance.x = p / solutionState.p;

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
            extent={{54,108},{-46,8}},
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
  end ExternalGas;

  model TerminalInflow "Molar pump of substance to system"

    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Chemical.Interfaces.ConditionalSubstanceFlow(useSubstanceFlowInput=false);

    import Chemical.Utilities.Types.SolutionChoice;

    parameter Chemical.Interfaces.Definition definition = dropOfCommons.DefaultSubstance
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Modelica.Units.SI.Time TC=dropOfCommons.TC "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    Chemical.Interfaces.Fore fore "Forwards port" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);
    outer Chemical.DropOfCommons dropOfCommons;

  initial equation
    u = u_start;

  equation
    fore.definition = definition;
    fore.solution_forwards = solutionState;
    connect(fore.solution_rearwards,inputSubstrateSolution);

    fore.n_flow = -q;

    TC * der(u) = fore.r;
    fore.state_forwards.u = u;
    fore.state_forwards.h = Chemical.Interfaces.Properties.molarEnthalpy(definition,solutionState);
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
  end TerminalInflow;

  model TerminalOutflow "Molar pump of substance from system"
    //  extends Chemical.Undirected.Boundaries.Internal.PartialTerminalRear;
    extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
    extends Chemical.Interfaces.ConditionalSubstanceFlow;

    import Chemical.Utilities.Types.SolutionChoice;

     Chemical.Interfaces.Rear rear "The substance"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

    parameter Modelica.Units.SI.Time TC=dropOfCommons.TC "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

  protected
    outer Chemical.DropOfCommons dropOfCommons "Chemical wide properties";
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_start;
  equation
    rear.n_flow = q;
    rear.solution_rearwards = solutionState;
    connect(rear.solution_forwards,inputSubstrateSolution);

    TC * der(u) = rear.r;
    rear.state_rearwards.u = u;
    rear.state_rearwards.h = Chemical.Interfaces.Properties.molarEnthalpy(rear.definition,solutionState);

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
  end TerminalOutflow;

  model Clearance "ClearanceFlow of whole solution"
    extends Internal.PartialTerminalRear;

    parameter Modelica.Units.SI.VolumeFlowRate Clearance
    "Physiological clearance of the substance";

  equation

    assert(Clearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

    rear.n_flow = substance.c * Clearance;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                          graphics={
          Rectangle(
            extent={{-100,-50},{100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Polygon(
            points={{-80,25},{80,0},{-80,-25},{-80,25}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Text(
            extent={{-150,-20},{150,20}},
            textString="%name",
            lineColor={128,0,255},
            origin={2,-74},
            rotation=180)}),
      Documentation(revisions="<html>
<p><i>2009-2018 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<h4><font color=\"#008000\">Bidirectional mass flow by concentration</font></h4>
<p>Possible field values: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0.1\"><tr>
<td></td>
<td><h4>forward flow</h4></td>
<td><h4>backward flow</h4></td>
</tr>
<tr>
<td><h4>solutionFlow</h4></td>
<td><p align=\"center\">&gt;=0</p></td>
<td><p align=\"center\">&lt;=0</p></td>
</tr>
<tr>
<td><h4>q_in.q</h4></td>
<td><p align=\"center\">=solutionFlow*q_in.conc</p></td>
<td><p align=\"center\">=-q_out.q</p></td>
</tr>
<tr>
<td><h4>q_out.q</h4></td>
<td><p align=\"center\">=-q_in.q</p></td>
<td><p align=\"center\">=solutionFlow*q_out.conc</p></td>
</tr>
</table>
<br/>
</html>"));
  end Clearance;

  model Degradation "Flow of whole solution"
    extends Internal.PartialTerminalRear;

    parameter Modelica.Units.SI.Time HalfTime
    "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

    Modelica.Units.SI.AmountOfSubstance n;

  equation

    n = substance.x*solutionState.n;
    rear.n_flow = n * (Modelica.Math.log(2)/HalfTime);

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                          graphics={
          Rectangle(
            extent={{-100,-50},{100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Polygon(
            points={{-80,25},{80,0},{-80,-25},{-80,25}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Text(
            extent={{-150,-20},{150,20}},
            textString="%name",
            lineColor={128,0,255},
            origin={2,-74},
            rotation=180)}),
      Documentation(revisions="<html>
<p><i>2009-2018 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<h4><font color=\"#008000\">Bidirectional mass flow by concentration</font></h4>
<p>Possible field values: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0.1\"><tr>
<td></td>
<td><h4>forward flow</h4></td>
<td><h4>backward flow</h4></td>
</tr>
<tr>
<td><h4>solutionFlow</h4></td>
<td><p align=\"center\">&gt;=0</p></td>
<td><p align=\"center\">&lt;=0</p></td>
</tr>
<tr>
<td><h4>q_in.q</h4></td>
<td><p align=\"center\">=solutionFlow*q_in.conc</p></td>
<td><p align=\"center\">=-q_out.q</p></td>
</tr>
<tr>
<td><h4>q_out.q</h4></td>
<td><p align=\"center\">=-q_in.q</p></td>
<td><p align=\"center\">=solutionFlow*q_out.conc</p></td>
</tr>
</table>
<br/>
</html>"));
  end Degradation;

  model BoundaryRear "Generic Boundary model (may act as source or sink)"



    parameter Chemical.Interfaces.Definition substanceDefinition=Chemical.Substances.Liquid.Unknown "Definition of the substance"
      annotation (choicesAllMatching=true, Dialog(group= "Substance"));

    parameter Boolean useSolution = false "Use input connector for solution?"
      annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Chemical solution"));
    parameter Chemical.Interfaces.SolutionState solutionParam = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible)
      annotation (Dialog(enable=not useSolution,group = "Chemical solution"));

    parameter Boolean usePotential = false "Use input connector for chemical potential"
      annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Substance"));
    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value"
      annotation (Dialog(enable=not usePotential, group= "Substance"));

    parameter Boolean useEnthalpy = false "Use input connector for molar enthalpy"
      annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Substance"));
    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
      annotation (Dialog(enable=not useEnthalpy, group= "Substance"));

    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation (Dialog(tab="Advanced"));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol")
                                                           if usePotential "Chemical potential input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,40},{0,80}}), iconTransformation(extent={{-40,40},{0,80}})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit="J/mol")  if  useEnthalpy "Enthalpy input connector"
      annotation (Placement(transformation(extent={{-40,-40},{0,0}}), iconTransformation(extent={{-40,-20},{0,20}})));
    Interfaces.Fore fore
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential u_rearwards=fore.state_rearwards.u;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

    Modelica.Units.SI.ChemicalPotential r;
    Chemical.Interfaces.SolutionState solutionState "State of chemical solution";

  public
    Chemical.Interfaces.SolutionPort solution(T=solutionState.T,p=solutionState.p,v=solutionState.v,n=solutionState.n,m=solutionState.m,V=solutionState.V,G=solutionState.G,Q=solutionState.Q,I=solutionState.I, i=0, dH=0, dV=0, nj=0, mj=0, Vj=0, Gj=0, Qj=0, Ij=0) if useSolution
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}}), iconTransformation(extent={{-30,-70},{-10,-50}})));
  equation

    connect(u0_var, u0);
    if not usePotential then
      u0 = u0_par;
    end if;

    connect(h0_var, h0);
    if not useEnthalpy then
       h0 = h0_par;
    end if;

    if not useSolution then
      solutionState=solutionParam;
    end if;

    der(fore.n_flow)*L = fore.r-r;

    //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
    r =.Chemical.Utilities.Internal.regStep(
        fore.n_flow,
        u0 - u_rearwards,
        0,
        n_flow_reg);

    fore.state_forwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);
    fore.solution_forwards=solutionState;
    fore.definition=substanceDefinition;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{0,76},{64,-84}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{0,80},{60,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{60,0},{84,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{60,80},{60,-80}},
            color={158,66,200},
            thickness=0.5),
          Line(points={{44,80},{44,-80}}, color={255,255,255}),
          Line(
            points={{28,80},{28,-80}},
            color={255,255,255},
            thickness=0.5),
          Line(
            points={{12,80},{12,-80}},
            color={255,255,255},
            thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>A undirected boundary that can act as source and sink, depending on the rest of the system. The Boundary_rear has to be connected to the rear end of your model and therefore has a fore port.</u>
<u>At positive massflow the fore port acts as an outlet and therefore the boundary_rear is a source.</u>
</html>"));
  end BoundaryRear;

  model BoundaryFore "Generic Boundary model (may act as source or sink)"

    parameter Boolean usePotential = false "Use input connector for potential?"
        annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Substance"));
    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value"
        annotation (Dialog(enable=not usePotential, group= "Substance"));

    parameter Boolean useEnthalpy = false "Use input connector for molar enthalpy"
        annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Substance"));
    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
        annotation (Dialog(enable=not useEnthalpy, group= "Substance"));


    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

    parameter Boolean useSolution = false "Use input connector for solution?"
      annotation ( Evaluate=true, HideResult=true, choices(checkBox=true), Dialog(group= "Chemical solution"));
    parameter Chemical.Interfaces.SolutionState solutionParam = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible)
      annotation (Dialog(enable=not useSolution,group = "Chemical solution"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol") if usePotential "ChemicalPotential input connector [Pa]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60}),
        iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/mol") if useEnthalpy "Enthalpy input connector"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,-20}),
        iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,0})));
    Chemical.Interfaces.Rear rear
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-80,-20},{-120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential u_forwards=rear.state_forwards.u;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

    Modelica.Units.SI.ChemicalPotential r;
    Chemical.Interfaces.SolutionState solutionState "State of chemical solution";

  public
    Chemical.Interfaces.SolutionPort solution(T=solutionState.T,p=solutionState.p,v=solutionState.v,n=solutionState.n,m=solutionState.m,V=solutionState.V,G=solutionState.G,Q=solutionState.Q,I=solutionState.I, i=0, dH=0, dV=0, nj=0, mj=0, Vj=0, Gj=0, Qj=0, Ij=0) if useSolution
      annotation (Placement(transformation(extent={{10,-70},{30,-50}}),   iconTransformation(extent={{10,-70},{30,-50}})));

  equation

    connect(u0_var, u0);
    if not usePotential then
      u0 = u0_par;
    end if;

    connect(h0_var, h0);
    if not useEnthalpy then
       h0 = h0_par;
    end if;

    if not useSolution then
      solutionState=solutionParam;
    end if;

    der(rear.n_flow)*L = rear.r-r;

    //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
    r = Chemical.Utilities.Internal.regStep(
        rear.n_flow,
        u0 - u_forwards,
        0,
        n_flow_reg);

    rear.state_rearwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);

    rear.solution_rearwards=solutionState;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{4,76},{-60,-84}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{0,80},{-60,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-60,0},{-84,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{-60,80},{-60,-80}},
            color={158,66,200},
            thickness=0.5),
          Line(points={{-44,80},{-44,-80}}, color={255,255,255}),
          Line(
            points={{-26,80},{-26,-80}},
            color={255,255,255},
            thickness=0.5),
          Line(
            points={{-12,80},{-12,-80}},
            color={255,255,255},
            thickness=1)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>A undirected boundary that can act as source and sink, depending on the rest of the system. The Boundary_fore has to be connected to the fore end of your model and therefore has a rear port.</u>
<u>At positive massflow the rear port acts as an inlet and therefore the boundary_fore is a sink.</u>
</html>"));
  end BoundaryFore;

  package Tests "Tests for the boundaries package"
    extends Modelica.Icons.ExamplesPackage;

    model TestSubstance
       extends Modelica.Icons.Example;
      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

      BoundaryRear boundaryRear(substanceDefinition=Chemical.Substances.Liquid.H2O, useSolution=false)
        annotation (Placement(transformation(extent={{-76,30},{-56,50}})));
      Substance substance(
        useRear=true,
        useFore=false,

        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{26,30},{46,50}})));
      Substance substance2(
        useRear=false,
        useFore=true,
        useSolution=true,
        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
      BoundaryFore boundaryFore annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      Substance substance1(
        useRear=false,
        useFore=true,
        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
      BoundaryRear boundaryRear1(substanceDefinition=Chemical.Substances.Liquid.H2O, useSolution=false)
        annotation (Placement(transformation(extent={{-74,50},{-54,70}})));
      Substance substance3(
        useRear=true,
        useFore=true,
        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{-24,50},{-4,70}})));
      BoundaryFore boundaryFore2
                                annotation (Placement(transformation(extent={{34,50},{54,70}})));
      BoundaryRear boundaryRear2(substanceDefinition=Chemical.Substances.Liquid.H2O, useSolution=false)
        annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      Substance substance4(
        useRear=true,
        useFore=false,

        useSolution=true,
        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
      BoundaryRear boundaryRear3(substanceDefinition=Chemical.Substances.Liquid.H2O, useSolution=false)
        annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      Substance substance5(
        useRear=true,
        useFore=true,
        useSolution=true,
        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
      BoundaryFore boundaryFore3
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryFore boundaryFore1 annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={44,82})));
      Substance solvent(
        substanceDefinition=Chemical.Substances.Liquid.H2O,
                        useFore=false, useSolution=true) annotation (Placement(transformation(extent={{66,-78},{86,-58}})));
      Substance substance6(substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{-96,70},{-76,90}})));
      BoundaryRear boundaryRear4(substanceDefinition=Chemical.Substances.Liquid.H2O, useSolution=false)
        annotation (Placement(transformation(extent={{-76,6},{-56,26}})));
      Substance substance7(
        useRear=true,
        useFore=false,
        useRearSolution=false,

        substanceDefinition=Chemical.Substances.Liquid.H2O) annotation (Placement(transformation(extent={{26,6},{46,26}})));
    equation
      connect(boundaryRear.fore, substance.rear) annotation (Line(
          points={{-56,40},{26,40}},
          color={158,66,200},
          thickness=0.5));
      connect(substance2.fore,boundaryFore. rear) annotation (Line(
          points={{-48,-16},{36,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(substance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}},        color={127,127,0}));
      connect(boundaryRear1.fore, substance3.rear) annotation (Line(
          points={{-54,60},{-24,60}},
          color={158,66,200},
          thickness=0.5));
      connect(substance3.fore, boundaryFore2.rear) annotation (Line(
          points={{-4,60},{34,60}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear2.fore, substance4.rear) annotation (Line(
          points={{-46,-72},{34,-72}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear3.fore, substance5.rear) annotation (Line(
          points={{-48,-44},{-18,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(substance5.fore, boundaryFore3.rear) annotation (Line(
          points={{2,-44},{40,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(substance5.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(substance4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(substance1.fore, boundaryFore1.rear) annotation (Line(
          points={{-52,82},{34,82}},
          color={158,66,200},
          thickness=0.5));
      connect(solution.solution, solvent.solution) annotation (Line(points={{60,-98.94},{60,-104},{70,-104},{70,-78}}, color={127,127,0}));
      connect(boundaryRear4.fore, substance7.rear) annotation (Line(
          points={{-56,16},{26,16}},
          color={158,66,200},
          thickness=0.5));
       annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=5,
          Interval=0.01,
          Tolerance=1e-06,
          __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
    end TestSubstance;

    model TestExternalSubstance
       extends Modelica.Icons.Example;
      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,6}})));



      ExternalSubstance
                externalSubstance1(
        useRear=true,
        useFore=false,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        FixedValue=1)     annotation (Placement(transformation(extent={{24,14},{44,34}})));

      ExternalSubstance externalSubstance2(
        useRear=false,
        useFore=true,
        useSolution=true,
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=10) annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));

      BoundaryFore boundaryFore annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      ExternalSubstance externalIdealGas(
        useRear=false,
        useFore=true,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        FixedValue=15) annotation (Placement(transformation(extent={{-70,72},{-50,92}})));
      BoundaryFore boundaryFore1
                                annotation (Placement(transformation(extent={{32,72},{52,92}})));
      BoundaryRear boundaryRear1(
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        useSolution=false) annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
      ExternalSubstance
                externalSubstance(
        useRear=true,
        useFore=true,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        FixedValue=10)    annotation (Placement(transformation(extent={{-26,42},{-6,62}})));

      BoundaryFore boundaryFore2
                                annotation (Placement(transformation(extent={{30,42},{50,62}})));
      BoundaryRear boundaryRear2(
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        useSolution=false) annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      ExternalSubstance
                externalSubstance4(
        useRear=true,        useFore=false,
        useSolution=true,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=1)     annotation (Placement(transformation(extent={{34,-82},{54,-62}})));

      BoundaryRear boundaryRear3(
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        useSolution=false) annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      ExternalSubstance
                externalSubstance3(
        useRear=true,
        useFore=true,
        useSolution=true,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=1)     annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));

      BoundaryFore boundaryFore3
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryRear boundaryRear4(
        substanceDefinition=Chemical.Substances.Liquid.Ethanol,
        useSolution=false) annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      Substance solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{70,-82},{90,-62}})));
    equation
      connect(externalSubstance2.fore, boundaryFore.rear) annotation (Line(
          points={{-48,-16},{36,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(externalSubstance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(externalIdealGas.fore, boundaryFore1.rear) annotation (Line(
          points={{-50,82},{32,82}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear1.fore, externalSubstance.rear) annotation (Line(
          points={{-56,52},{-26,52}},
          color={158,66,200},
          thickness=0.5));
      connect(externalSubstance.fore, boundaryFore2.rear) annotation (Line(
          points={{-6,52},{30,52}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear2.fore, externalSubstance4.rear) annotation (Line(
          points={{-46,-72},{34,-72}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear3.fore, externalSubstance3.rear) annotation (Line(
          points={{-48,-44},{-18,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(externalSubstance3.fore, boundaryFore3.rear) annotation (Line(
          points={{2,-44},{40,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(externalSubstance3.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(externalSubstance4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(boundaryRear4.fore, externalSubstance1.rear) annotation (Line(
          points={{-54,24},{24,24}},
          color={158,66,200},
          thickness=0.5));
      connect(solution.solution, solvent.solution) annotation (Line(points={{60,-98.94},{60,-104},{74,-104},{74,-82}}, color={127,127,0}));
       annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
    end TestExternalSubstance;

    model TestExternalGas
      import Chemical;
       extends Modelica.Icons.Example;
      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

      ExternalGas
                externalGas1(
        useRear=true,
        useFore=false,
        PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{24,14},{44,34}})));
      ExternalGas externalGas2(
        useRear=false,
        useFore=true,
        useSolution=true,
        substanceDefinition=Chemical.Substances.Gas.H2O,
        PartialPressure(displayUnit="mmHg") = 1333.22387415) annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
      BoundaryFore boundaryFore annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      ExternalGas externalIdealGas(
        useRear=false,
        useFore=true,
        substanceDefinition=Chemical.Substances.Gas.H2O,
        PartialPressure(displayUnit="mmHg") = 1999.835811225) annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
      BoundaryFore boundaryFore1
                                annotation (Placement(transformation(extent={{32,72},{52,92}})));
      BoundaryRear boundaryRear1(
        substanceDefinition=Chemical.Substances.Gas.H2O,
        useSolution=false) annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
      ExternalGas
                externalGas(
        useRear=true,
        useFore=true,
        PartialPressure(displayUnit="mmHg") = 1333.22387415)
                          annotation (Placement(transformation(extent={{-24,42},{-4,62}})));
      BoundaryFore boundaryFore2
                                annotation (Placement(transformation(extent={{30,42},{50,62}})));
      BoundaryRear boundaryRear2(
        substanceDefinition=Chemical.Substances.Gas.H2O,
        useSolution=false) annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      ExternalGas
                externalGas4(
        useRear=true,        useFore=false,
        useSolution=true,
        PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
      BoundaryRear boundaryRear3(
        substanceDefinition=Chemical.Substances.Gas.H2O,
        useSolution=false) annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      ExternalGas
                externalGas3(
        useRear=true,
        useFore=true,
        useSolution=true,
        PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
      BoundaryFore boundaryFore3
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryRear boundaryRear4(
        substanceDefinition=Chemical.Substances.Gas.H2O,
        useSolution=false) annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      Chemical.Boundaries.Substance solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{66,-76},{86,-56}})));
    equation
      connect(externalGas2.fore, boundaryFore.rear) annotation (Line(
          points={{-48,-16},{36,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(externalGas2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(externalIdealGas.fore, boundaryFore1.rear) annotation (Line(
          points={{-52,82},{32,82}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear1.fore, externalGas.rear) annotation (Line(
          points={{-56,52},{-24,52}},
          color={158,66,200},
          thickness=0.5));
      connect(externalGas.fore, boundaryFore2.rear) annotation (Line(
          points={{-4,52},{30,52}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear2.fore, externalGas4.rear) annotation (Line(
          points={{-46,-72},{34,-72}},
          color={158,66,200},
          thickness=0.5));
      connect(boundaryRear3.fore, externalGas3.rear) annotation (Line(
          points={{-48,-44},{-18,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(externalGas3.fore, boundaryFore3.rear) annotation (Line(
          points={{2,-44},{40,-44}},
          color={158,66,200},
          thickness=0.5));
      connect(externalGas3.solution, solution.solution) annotation (Line(points={{-14,-54},{-14,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(externalGas4.solution, solution.solution) annotation (Line(points={{38,-82},{38,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(boundaryRear4.fore, externalGas1.rear) annotation (Line(
          points={{-54,24},{24,24}},
          color={158,66,200},
          thickness=0.5));
      connect(solution.solution, solvent.solution) annotation (Line(points={{60,-98.94},{60,-104},{70,-104},{70,-76}}, color={127,127,0}));
       annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
    end TestExternalGas;

    model TestBoundaries "Tests for the rear and fore boundary"
      extends Modelica.Icons.Example;

      BoundaryRear boundary_rear(
        u0_par=100000,
        fore(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-28,82})));
      BoundaryFore boundary_fore(
        usePotential=true,
        u0_par=110000) annotation (Placement(transformation(extent={{22,72},{42,92}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-88,72},{-68,92}})));
      Modelica.Blocks.Sources.Step step(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{62,76},{50,88}})));
      BoundaryFore boundary_fore1(
        usePotential=true,
        u0_par=110000) annotation (Placement(transformation(extent={{22,46},{42,66}})));
      Modelica.Blocks.Sources.Step step1(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{62,50},{50,62}})));
      BoundaryRear boundary_rear1(
        usePotential=true,
        u0_par=100000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-26,30})));
      Modelica.Blocks.Sources.Step step2(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{-52,24},{-40,36}})));
      BoundaryRear boundary_rear2(
        useSolution=true,
        u0_par=100000,
        fore(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-22})));
      BoundaryFore boundary_fore2(usePotential=true, u0_par=110000)
                       annotation (Placement(transformation(extent={{-10,-32},{10,-12}})));
      Modelica.Blocks.Sources.Step step3(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{30,-28},{18,-16}})));
      BoundaryFore boundary_fore3(usePotential=true, u0_par=110000)
                       annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
      Modelica.Blocks.Sources.Step step4(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{64,-56},{52,-44}})));
      BoundaryRear boundary_rear3(
        substanceDefinition=Chemical.Substances.Liquid.H2O,
        useSolution=true,
        usePotential=true,
        u0_par=100000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={26,-78})));
      Modelica.Blocks.Sources.Step step5(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{-10,-90},{2,-78}})));
      Solution solution annotation (Placement(transformation(extent={{-98,-98},{102,0}})));
      TerminalOutflow terminalOutflow annotation (Placement(transformation(extent={{32,20},{52,40}})));
      TerminalInflow terminalInflow                    annotation (Placement(transformation(extent={{-48,50},{-28,70}})));
      TerminalOutflow terminalOutflow1(SubstanceFlow=0) annotation (Placement(transformation(extent={{66,-86},{86,-66}})));
      TerminalInflow terminalInflow1(solutionFrom=Chemical.Utilities.Types.SolutionChoice.SolutionPort)
                                     annotation (Placement(transformation(extent={{-48,-54},{-28,-34}})));
    equation
      connect(step.y, boundary_fore.u0_var)
        annotation (Line(points={{49.4,82},{42,82},{42,88},{34,88}},
                                                       color={0,0,127}));
      connect(boundary_fore.rear, boundary_rear.fore) annotation (Line(
          points={{22,82},{-18,82}},
          color={158,66,200},
          thickness=0.5));
      connect(step1.y, boundary_fore1.u0_var) annotation (Line(points={{49.4,56},{42,56},{42,62},{34,62}},
                                                                                             color={0,0,127}));
      connect(step2.y,boundary_rear1.u0_var)  annotation (Line(points={{-39.4,30},{-34,30},{-34,24},{-28,24}},
                                                                                               color={0,0,127}));
      connect(step3.y, boundary_fore2.u0_var) annotation (Line(points={{17.4,-22},{10,-22},{10,-16},{2,-16}}, color={0,0,127}));
      connect(boundary_fore2.rear, boundary_rear2.fore) annotation (Line(
          points={{-10,-22},{-50,-22}},
          color={158,66,200},
          thickness=0.5));
      connect(step4.y,boundary_fore3. u0_var) annotation (Line(points={{51.4,-50},{44,-50},{44,-44},{36,-44}},
                                                                                             color={0,0,127}));
      connect(step5.y,boundary_rear3.u0_var)  annotation (Line(points={{2.6,-84},{24,-84}},    color={0,0,127}));
      connect(boundary_rear2.solution, solution.solution) annotation (Line(points={{-62,-16},{-62,-102},{62,-102},{62,-97.02}}, color={127,127,0}));
      connect(boundary_rear3.solution,solution. solution) annotation (Line(points={{24,-72},{24,-102},{62,-102},{62,-97.02}},
                                                                                                                  color={127,127,0}));
      connect(terminalOutflow.rear, boundary_rear1.fore) annotation (Line(
          points={{32,30},{-16,30}},
          color={158,66,200},
          thickness=0.5));
      connect(boundary_rear3.fore, terminalOutflow1.rear)
        annotation (Line(
          points={{36,-78},{36,-76},{66,-76}},
          color={158,66,200},
          thickness=0.5));
      connect(terminalInflow.fore, boundary_fore1.rear)
        annotation (Line(
          points={{-28,60},{12,60},{12,56},{22,56}},
          color={158,66,200},
          thickness=0.5));
      connect(terminalInflow1.fore, boundary_fore3.rear)
        annotation (Line(
          points={{-28,-44},{-28,-46},{24,-46},{24,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(terminalInflow1.solution, solution.solution) annotation (Line(points={{-44,-54},{-44,-97.02},{62,-97.02}}, color={127,127,0}));
      annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
    end TestBoundaries;

    annotation (Documentation(info="<html>
<u>Tests for the boundaries package.</u>
</html>"));
  end Tests;

  package Internal "Partials and Internal functions"
  extends Modelica.Icons.InternalPackage;

    partial model PartialBoundaryBase "Base boundary"


      outer Modelica.Fluid.System system "System wide properties";


      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L
       annotation(HideResult=true, Dialog(tab = "Advanced"));

      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation(HideResult=true, Dialog(tab="Advanced"));

     protected


      Modelica.Units.SI.MolarFlowRate n_flow "Molar change of the amount of base substance";
      Modelica.Units.SI.EnthalpyFlowRate h_flow "Change of enthalpy";

      outer Chemical.DropOfCommons dropOfCommons "Chemical wide properties";


       //if port.n_flow > 0 -> it is sink (r=medium.u-u_in) else it is source (r=0)
      Modelica.Units.SI.ChemicalPotential r_rear_intern=Chemical.Utilities.Internal.regStep(
                n_flow_rear,
                state_out.u - state_in_rear.u,
                0,
                n_flow_reg);
      Modelica.Units.SI.ChemicalPotential r_fore_intern=Chemical.Utilities.Internal.regStep(
                n_flow_fore,
                state_out.u - state_in_fore.u,
                0,
                n_flow_reg);
      // dont regstep variables that are only in der(state), to increase accuracy
      Modelica.Units.SI.EnthalpyFlowRate h_flow_rear=(if n_flow_rear >= 0 then state_in_rear.h else state_out.h)*n_flow_rear;
      Modelica.Units.SI.EnthalpyFlowRate h_flow_fore=(if n_flow_fore >= 0 then state_in_fore.h else state_out.h)*n_flow_fore;

      Modelica.Units.SI.ChemicalPotential r_rear_port;
      Modelica.Units.SI.ChemicalPotential r_fore_port;
      Modelica.Units.SI.MolarFlowRate n_flow_rear;
      Modelica.Units.SI.MolarFlowRate n_flow_fore;

      Chemical.Interfaces.SubstanceStateInput state_in_rear;
      Chemical.Interfaces.SubstanceStateInput state_in_fore;
      Chemical.Interfaces.SubstanceState state_out;

    equation

      der(n_flow_rear)*L = r_rear_port - r_rear_intern;
      der(n_flow_fore)*L = r_fore_port - r_fore_intern;


      n_flow = n_flow_rear + n_flow_fore;
      h_flow = h_flow_rear + h_flow_fore;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialBoundaryBase;

    partial model PartialBoundary "Boundary enabling rear and fore connector"
      extends PartialBoundaryBase;

      parameter Chemical.Interfaces.Definition substanceDefinition=dropOfCommons.DefaultSubstance "Definition of the substance"
        annotation (choicesAllMatching=true, Dialog(enable=not useRear, group = "Substance"));

      parameter Boolean useRear = false "Use rearwards conector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));
      parameter Boolean useFore = true "Use forwards connector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));


      Chemical.Interfaces.Rear rear(
        n_flow=n_flow_rear,
        r=r_rear_port,
        state_rearwards=state_out,
        solution_rearwards=solutionState) if useRear annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

      Chemical.Interfaces.Fore fore(
        n_flow=n_flow_fore,
        r=r_fore_port,
        state_forwards=state_out,
        solution_forwards=solutionState,
        definition=substanceDefinitionVar) if useFore annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));


    protected



      Chemical.Interfaces.Definition substanceDefinitionVar;         //substanceDefinition;
      Chemical.Interfaces.DefinitionInput definitionVar;//substanceDefinition;
      Chemical.Interfaces.SolutionState solutionState;


    equation
      substanceDefinitionVar = definitionVar;


      connect(state_in_rear,rear.state_forwards);
      connect(state_in_fore,fore.state_rearwards);


      if not useRear then
        r_rear_port = 0;
        n_flow_rear = 0;
        state_in_rear.h = 0;
        substanceDefinitionVar = substanceDefinition;
      else
        connect(rear.definition, definitionVar);
      end if;

      if not useFore then
        r_fore_port = 0;
        n_flow_fore = 0;
        state_in_fore.h = 0;
      end if;


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialBoundary;

    partial model PartialSubstance "Chemical substance base"
      extends PartialBoundary;


      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance"
       annotation (Dialog(group="Substance"));

      Modelica.Units.SI.AmountOfSubstance n
        "Amount of base molecules inside all clusters in compartment";

      Chemical.Interfaces.Properties.SubstanceProperties substance(
        definition=substanceDefinitionVar,
        solutionState=solutionState,
        FixedDefinition=not useRear,
        definitionParam=substanceDefinition,
        amountOfBaseMolecules=n,
        m_start=m_start,
        n_flow=n_flow,
        h_flow=h_flow);

      parameter Boolean useSolution = false "Use solution port?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Chemical solution"));
      parameter Boolean useRearSolution = useRear and (not useSolution)  "Use solution from rear?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(enable=not useSolution and useRear, group="Chemical solution"));

      parameter Chemical.Interfaces.SolutionState solutionParam = Chemical.Interfaces.SolutionState(phase=Chemical.Interfaces.Phase.Incompressible) "Constant chemical solution state if not from rear or input"
          annotation (HideResults=useSolution or (useRearSolution and useRear), Dialog(enable=not useSolution and (not useRear or not useRearSolution), group="Chemical solution"));
      Chemical.Interfaces.SolutionPort solution(
          T=solutionPortState.T,
          p=solutionPortState.p,
          v=solutionPortState.v,
          n=solutionPortState.n,
          m=solutionPortState.m,
          V=solutionPortState.V,
          G=solutionPortState.G,
          Q=solutionPortState.Q,
          I=solutionPortState.I,
          i=substance.i,
          dH=substance.dH,
          dV=substance.dV,
          nj=substance.nj,
          mj=substance.mj,
          Vj=substance.Vj,
          Gj=substance.Gj,
          Qj=substance.Qj,
          Ij=substance.Ij)
            if useSolution "To connect substance with solution, where is pressented"
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));
        Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
            "Concentration";
        Modelica.Units.SI.MassConcentration M(displayUnit="g/l")
            "Mass concentration";
        Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
            "Molality";
        Modelica.Units.SI.MoleFraction x "Mole fraction";
        Modelica.Units.SI.MassFraction X "Mass fraction";
    protected
         outer Chemical.DropOfCommons dropOfCommons "Chemical wide properties";
         Chemical.Interfaces.SolutionState solutionPortState;

         Chemical.Interfaces.SolutionStateInput inputSubstrateSolution=solutionState if (useRearSolution and useRear and not useSolution);
    equation



      c = substance.c;
      b = substance.b;
      M = substance.M;
      x = substance.x;
      X = substance.X;
      state_out.u = substance.u;
      state_out.h = substance.h;

      connect(rear.solution_forwards,inputSubstrateSolution);

      if (useSolution and not useRearSolution) or (not useSolution) then
        solutionState=solutionPortState;
      end if;
      if not useSolution and not useRearSolution then
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
        annotation (choicesAllMatching=true, HideResults=useSolution or useRear, Dialog(enable=not useSolution and not useRear),
                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));


    end PartialSubstance;

    partial model PartialTerminalRear
     extends Chemical.Interfaces.PartialSolutionSensor(solutionFrom = SolutionChoice.Parameter);
     import Chemical.Utilities.Types.SolutionChoice;

     outer Modelica.Fluid.System system "System wide properties";


      Chemical.Interfaces.Rear rear "The substance"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

      Chemical.Interfaces.Properties.SubstanceProperties substance(
        SolutionObserverOnly=true,
        definitionParam=Chemical.Substances.Liquid.H2O,
        definition(
          data=rear.definition.data,
          SelfClustering=false,
          SelfClustering_dH=0,
          SelfClustering_dS=0),
        solutionState=solutionState,
        FixedDefinition=false,
        m_start=1,
        n_flow=0,
        h_flow=0);

      parameter Modelica.Units.SI.Time TC=dropOfCommons.TC "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    protected
      Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

      outer Chemical.DropOfCommons dropOfCommons "Chemical wide properties";
    initial equation
      u = u_start;
    equation
      rear.solution_rearwards = solutionState;
      connect(rear.solution_forwards,inputSubstrateSolution);

     substance.u = rear.state_forwards.u;

      TC * der(u) = rear.r;
      rear.state_rearwards.u = u;
      rear.state_rearwards.h = Chemical.Interfaces.Properties.molarEnthalpy(rear.definition,solutionState);

     annotation (
       Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialTerminalRear;

    partial model ConditionalSolutionFlow "Input of solution molar flow vs. parametric solution molar flow"

      parameter Boolean useSolutionFlowInput = false
      "=true, if solution flow is provided via input"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
              Dialog(group="Conditional inputs", __Dymola_compact=true));

    parameter Modelica.Units.SI.VolumeFlowRate SolutionFlow=0
      "Volume flow rate of the solution if useSolutionFlowInput=false"
      annotation (HideResult=true, Dialog(enable=not useSolutionFlowInput));

      Modelica.Blocks.Interfaces.RealInput solutionFlow(start=SolutionFlow, final unit="m3/s")=
         volumeFlow if useSolutionFlowInput
         annotation ( HideResult=true, Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,40})));

    Modelica.Units.SI.VolumeFlowRate volumeFlow "Current solution volume flow";

    equation
      if not useSolutionFlowInput then
        volumeFlow = SolutionFlow;
      end if;

    end ConditionalSolutionFlow;

    package Types
      type ConcentrationQuantities = enumeration(
          c_molpm3 "Concentration (mmol/L)",
          X_kgpkg "Mass fraction (kg/kg)",
          b_molpkg "Molality (mol/kg)",
          x_molpmol "Mole fraction (mol/mol)",
          p_Pa "Partial pressure (Pa)",
          p_kPa "Partial pressure (kPa)",
          p_mmHg "Partial pressure (mmHg)",
          p_bar "Partial pressure (bar)");
    end Types;

    function getUnit "Returns unit of input quantity"
      extends Modelica.Icons.Function;

      input Types.ConcentrationQuantities quantity;
      output String unit;

    algorithm

      if quantity == Types.ConcentrationQuantities.c_molpm3 then
        unit := "mol/m3";
      elseif quantity == Types.ConcentrationQuantities.X_kgpkg then
        unit := "kg/kg";
      elseif quantity == Types.ConcentrationQuantities.b_molpkg then
        unit := "mol/kg";
      elseif quantity == Types.ConcentrationQuantities.x_molpmol then
        unit := "mol/mol";
      elseif quantity == Types.ConcentrationQuantities.p_Pa then
        unit := "Pa";
      elseif quantity == Types.ConcentrationQuantities.p_kPa then
        unit := "kPa";
      elseif quantity == Types.ConcentrationQuantities.p_mmHg then
        unit := "mmHg";
      elseif quantity == Types.ConcentrationQuantities.p_bar then
        unit := "bar";
      else
        unit :="";
      end if;

      annotation (Documentation(info="<html>
<p>Helper function to get the unit for a quantity.</p>
</html>"));
    end getUnit;
  annotation (Documentation(info="<html>
<p>This package contains all internal functions, partials, and other (e.g. experimental) models for the Boundaries package.</p>
</html>"));
  end Internal;

annotation (Documentation(info="<html>
<u>
Boundary models for undirected chemical simulation.
</u>
</html>"));
end Boundaries;
