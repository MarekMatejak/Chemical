within Chemical;
package Boundaries "Boundary models for undirected chemical simulation"
  extends Modelica.Icons.SourcesPackage;

  model Substance "Substance in solution"
    extends Icons.Substance;
    extends Internal.PartialSubstance(
      useFore=false,
      useSolution=false,
      useRear=false,
      m_start=if use_mass_start then mass_start else
       amountOfSubstance_start*stateOfMatter.molarMassOfBaseMolecule(substanceData));

    import Chemical.Utilities.Types.InitializationUndirectedSubstance;

    parameter Boolean use_mass_start=true  "prefere state as mass, otherwise amountOfSubstance"
      annotation (Evaluate=true, choices(checkBox=true));

    parameter Modelica.Units.SI.Mass mass_start=1 "Initial mass of the substance"
      annotation (HideResult=not use_mass_start, Dialog(group="Initialization", enable=use_mass_start));

    parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=1
    "Initial amount of substance base molecules"
      annotation ( Dialog(group="Initialization", enable=(not use_mass_start)));

    Modelica.Units.SI.Mass mass=n*molarMassOfBaseMolecule "Mass";

    parameter InitializationUndirectedSubstance initAmount = Chemical.Utilities.Types.InitializationUndirectedSubstance.state "Initialization"
      annotation(Dialog(tab="Initialization"));

    parameter Modelica.Units.SI.MolarFlowRate change_start=1
    "Initial molar change of substance base molecules"
      annotation ( Dialog(tab="Initialization", enable=(initAmount == InitializationMethods.derivative)));

  protected

    Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceDataVar);

    Real logn(stateSelect=StateSelect.prefer, start=log(amountOfSubstance_start))
    "Natural logarithm of the amount of base molecules in solution";

    Real logm(stateSelect=StateSelect.prefer, start=log(mass_start))
    "Natural logarithm of the substance mass in solution";

  initial equation

    if initAmount == InitializationUndirectedSubstance.steadyStateForwards then
      substance.u = state_in_rear.u;
    elseif initAmount == InitializationUndirectedSubstance.steadyStateRearwards then
      substance.u = state_in_fore.u;
    elseif initAmount == InitializationUndirectedSubstance.state and not use_mass_start then
      logn=log(amountOfSubstance_start);
    elseif initAmount == InitializationUndirectedSubstance.state and use_mass_start then
      logm=log(mass_start);
    elseif initAmount == InitializationUndirectedSubstance.derivative then
      n_flow = change_start;
    end if;

  equation

    //The main accumulation equation is "der(n)=n_flow"
    // However, the numerical solvers can handle it during equilibration of chemical potential in form of log(n) much better. :-)
    if use_mass_start then
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

  model ElectronSource "Electron transfer from the solution to electric circuit"
    extends Icons.ElectronTransfer;

    Chemical.Interfaces.Fore fore(
      r=r_out,
      n_flow=n_flow,
      state_forwards(u=u, h=h),
      solution(
        T=solution.T,
        p=solution.p,
        v=solution.v,
        n=solution.n,
        m=solution.m,
        V=solution.V,
        Q=solution.Q,
        I=solution.I,
        G=solution.G),
      definition=substanceData) "The substance exiting" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
        Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
            extent={{-10,88},{10,108}})));

    Chemical.Interfaces.SolutionPort solution(
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

    parameter Chemical.Interfaces.Incompressible.SubstanceDataParameters substanceData=Chemical.Substances.Electrone_solid() "Definition of the substance";

    Real r_out, h;

    Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

    Modelica.Units.SI.ElectricPotential electricPotential "Electric potential of the solution";

    Modelica.Units.SI.Temperature temperature "Temperature of the solution";

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L;
    parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial electro-chemical potential";

    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Chemical.Interfaces.InputSubstanceState state_in_fore;

  initial equation
    u = u_0;
  equation

    connect(state_in_fore,fore.state_rearwards);

    //electric
    pin.v = electricPotential;
    pin.i + substanceData.z*Modelica.Constants.F*n_flow + solution.i = 0;

    /*
  These equations :
  
  u = Interfaces.Incompressible.chemicalPotentialPure(
    substanceData,
    temperature,
    pressure,
    electricPotential,
    moleFractionBasedIonicStrength) + (Modelica.Constants.R*temperature)*log(a) + substanceData.z*Modelica.Constants.F*electricPotential;
  uRT = u/(Modelica.Constants.R*temperature);
  h = Interfaces.Incompressible.molarEnthalpy(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
  
  ... are simplified as:
  */
    u = substanceData.z*Modelica.Constants.F*electricPotential;
    h = substanceData.z*Modelica.Constants.F*electricPotential;

    // Bounsaries.Source - u adaptation
    der(u)*L = r_out;

    //solution changes

    temperature = solution.T;
    electricPotential = solution.v;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Text(
            extent={{-146,-44},{154,-84}},
            textString="%name",
            lineColor={128,0,255})}),
      Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ElectronSource;

  model ElectronSink "Electron transfer to an electric circuit"
    extends Icons.ElectronTransfer;

    Chemical.Interfaces.Rear rear(n_flow=n_flow) "Chemical electron inlet"
      annotation (Placement(transformation(extent={{110,-10},{90,10}}), iconTransformation(extent={{110,-10},{90,10}})));

    Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
        Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
            extent={{-10,88},{10,108}})));

    Chemical.Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

   parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L;

    Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

    Modelica.Units.SI.ElectricPotential electricPotential "Electric potential of the solution";

    Modelica.Units.SI.Temperature temperature "Temperature of the solution";

    Modelica.Units.SI.ChemicalPotential u;

   parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation (Dialog(tab="Advanced"));

  protected

    parameter Chemical.Interfaces.Incompressible.SubstanceData substanceData = Chemical.Substances.Electrone_solid() "Definition of the substance";

    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential r_rear_intern=Chemical.Utilities.Internal.regStep(
              n_flow,
              u - rear.state_forwards.u,
              0,
              n_flow_reg);

  initial equation
    u = rear.state_forwards.u;
  equation

    rear.state_rearwards.u=u;
    rear.state_rearwards.h=u;

  //  substanceData=rear.definition;

    der(rear.n_flow)*L = rear.r - r_rear_intern;

    //electric
    pin.v = electricPotential;
    pin.i + substanceData.z*Modelica.Constants.F*n_flow + solution.i = 0;

    /*
  These equations :
  
  u = Interfaces.Incompressible.chemicalPotentialPure(
    substanceData,
    temperature,
    pressure,
    electricPotential,
    moleFractionBasedIonicStrength) + (Modelica.Constants.R*temperature)*log(a) + substanceData.z*Modelica.Constants.F*electricPotential;
  uRT = u/(Modelica.Constants.R*temperature);
  
  ... are simplified as:
  */
    u = substanceData.z*Modelica.Constants.F*electricPotential;

    //solution changes
    solution.dH = 0;
    solution.dV = 0;

    //extensive properties of the solution
    solution.nj=0;
    solution.mj=0;
    solution.Vj=0;
    solution.Gj=0;
    solution.Qj=0;
    solution.Ij=0;

    temperature = solution.T;
    electricPotential = solution.v;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Text(
            extent={{-146,-44},{154,-84}},
            textString="%name",
            lineColor={128,0,255})}),
      Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ElectronSink;

  model ExternalSubstance "Constant source of molar concentration"
    extends Internal.PartialSubstance(
      m_start=1,
      substance(SolutionObserverOnly=true));

   /* parameter stateOfMatter.SubstanceDataParameters substanceData
 "Definition of the substance"
    annotation (choicesAllMatching = true, Dialog(enable=not useRear));
*/
    parameter Chemical.Boundaries.Internal.Types.ConcentrationQuantities quantity "Concentration quantity";

    parameter Real FixedValue = 1e-8
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
      value = 1*((substance.x * solutionState.n)/solutionState.m)/stateOfMatter.specificAmountOfParticles(substanceDataVar,
     solutionState.T,
     solutionState.p,
     solutionState.v,
     solutionState.I);
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
      redeclare package stateOfMatter = gasModel "Ideal Gas from MSL",
      m_start=1,
      substance(SolutionObserverOnly=true));

     replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
      Chemical.Interfaces.StateOfMatter "Gas substance model"
      annotation (choices(
        choice(redeclare package gasModel =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package gasModel =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package gasModel =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));
  /*
   parameter gasModel.SubstanceDataParameters substanceData
 "Definition of the substance"
    annotation (choicesAllMatching = true, Dialog(enable=not useRear));
*/
    parameter Boolean usePartialPressureInput = false
    "=true, if fixed partial pressure is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Pressure PartialPressure=1e-05
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
    extends Chemical.Boundaries.Internal.PartialSolutionSensor         (useSolutionFromRear=false);
    extends Chemical.Interfaces.ConditionalSubstanceFlow(useSubstanceFlowInput=false);

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

  //  parameter Chemical.Interfaces.SolutionStateParameters solutionState;
    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    Chemical.Interfaces.Fore fore "Forwards port" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_start;

  equation
    fore.definition = substanceData;
    fore.solution = solutionState;

    fore.n_flow = -q;

    TC * der(u) = fore.r;
    fore.state_forwards.u = u;
    fore.state_forwards.h = stateOfMatter.molarEnthalpy(substanceData,solutionState.T);
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
    extends Chemical.Interfaces.ConditionalSubstanceFlow;

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

    Chemical.Interfaces.Rear rear(redeclare package stateOfMatter = stateOfMatter) "The substance"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_start;
  equation
    rear.n_flow = q;

    TC * der(u) = rear.r;
    rear.state_rearwards.u = u;
    rear.state_rearwards.h = rear.stateOfMatter.molarEnthalpy(rear.definition,rear.solution.T);

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

    n = substance.x*rear.solution.n;
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

    parameter Chemical.Interfaces.SolutionStateParameters solutionState
      annotation (Dialog(enable=not solutionFromInput));
    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Boolean solutionFromInput = false "Use input connector for solution?";
    parameter Boolean potentialFromInput = false "Use input connector for chemical potential";
    parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";

    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
      annotation (Dialog(enable=not enthalpyFromInput));
    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value" annotation (Dialog(enable=not potentialFromInput));

    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation (Dialog(tab="Advanced"));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol")
                                                           if potentialFromInput "Chemical potential input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,40},{0,80}}), iconTransformation(extent={{-40,40},{0,80}})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit="J/mol")  if  enthalpyFromInput "Enthalpy input connector"
      annotation (Placement(transformation(extent={{-40,-40},{0,0}}), iconTransformation(extent={{-40,-20},{0,20}})));
    Interfaces.Fore fore(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential u_rearwards=fore.state_rearwards.u;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

    Modelica.Units.SI.ChemicalPotential r;
    Chemical.Interfaces.SolutionState s "State of chemical solution";

  public
    Chemical.Interfaces.SolutionPort solution(T=s.T,p=s.p,v=s.v,n=s.n,m=s.m,V=s.V,G=s.G,Q=s.Q,I=s.I, i=0, dH=0, dV=0, nj=0, mj=0, Vj=0, Gj=0, Qj=0, Ij=0) if solutionFromInput
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}}), iconTransformation(extent={{-30,-70},{-10,-50}})));
  equation

    connect(u0_var, u0);
    if not potentialFromInput then
      u0 = u0_par;
    end if;

    connect(h0_var, h0);
    if not enthalpyFromInput then
       h0 = h0_par;
    end if;

    if not solutionFromInput then
      s=solutionState;
    end if;

    der(fore.n_flow)*L = fore.r-r;

    //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
    r =.Chemical.Utilities.Internal.regStep(
        fore.n_flow,
        u0 - u_rearwards,
        0,
        n_flow_reg);

    fore.state_forwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);
    fore.solution=s;
    fore.definition=substanceData;

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

    parameter Boolean potentialFromInput = false "Use input connector for potential?";
    parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";

    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "molar enthalpy set value"
      annotation (Dialog(enable=not enthalpyFromInput));
    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "ChemicalPotential set value" annotation (Dialog(enable=not potentialFromInput));

    parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
      annotation (Dialog(tab="Advanced"));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the boundary" annotation (Dialog(tab="Advanced"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol") if potentialFromInput "ChemicalPotential input connector [Pa]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60}),
        iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,60})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/mol") if enthalpyFromInput "Enthalpy input connector"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,-20}),
        iconTransformation(extent={{-20,-20},{20,20}}, rotation=180, origin={20,0})));
    Chemical.Interfaces.Rear rear(redeclare package stateOfMatter = stateOfMatter)
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-80,-20},{-120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential u_forwards=rear.state_forwards.u;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

    Modelica.Units.SI.ChemicalPotential r;

  equation

    connect(u0_var, u0);
    if not potentialFromInput then
      u0 = u0_par;
    end if;

    connect(h0_var, h0);
    if not enthalpyFromInput then
       h0 = h0_par;
    end if;

    der(rear.n_flow)*L = rear.r-r;

    //if port.n_flow > 0 -> it is sink (r=u_set-u_in) else it is source (r=0)
    r =.Chemical.Utilities.Internal.regStep(
        rear.n_flow,
        u0 - u_forwards,
        0,
        n_flow_reg);

    rear.state_rearwards = Chemical.Interfaces.SubstanceState(u=u0,h=h0);

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

      BoundaryRear boundaryRear(
        substanceData=Chemical.Substances.Water_liquid(),
        solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-76,14},{-56,34}})));
      Substance substance(
        useRear=true,
        useFore=false,

        substanceData=Chemical.Substances.Water_liquid())
                          annotation (Placement(transformation(extent={{24,12},{44,32}})));
      Substance substance2(
        useRear=false,
        useFore=true,
        useSolution=true,
        substanceData=Chemical.Substances.Water_liquid())
                                                annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
      BoundaryFore boundaryFore annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      Substance substance1(useRear=false,
        useFore=true,                     substanceData=Chemical.Substances.Water_liquid())
                                                annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
      BoundaryRear boundaryRear1(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-78,42},{-58,62}})));
      Substance substance3(
        useRear=true,
        useFore=true,      substanceData=Chemical.Substances.Water_liquid())
                          annotation (Placement(transformation(extent={{-28,42},{-8,62}})));
      BoundaryFore boundaryFore2
                                annotation (Placement(transformation(extent={{30,42},{50,62}})));
      BoundaryRear boundaryRear2(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      Substance substance4(
        useRear=true,
        useFore=false,

        useSolution=true,
        substanceData=Chemical.Substances.Water_liquid())
                          annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
      BoundaryRear boundaryRear3(substanceData=Chemical.Substances.Water_liquid(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      Substance substance5(
        useRear=true,
        useFore=true,
        useSolution=true,                          substanceData=Chemical.Substances.Water_liquid())
                          annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
      BoundaryFore boundaryFore3
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryFore boundaryFore1 annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={44,82})));
      Substance solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{66,-78},{86,-58}})));
    equation
      connect(boundaryRear.fore, substance.rear) annotation (Line(
          points={{-56,24},{-16,24},{-16,22},{24,22}},
          color={158,66,200},
          thickness=0.5));
      connect(substance2.fore,boundaryFore. rear) annotation (Line(
          points={{-48,-16},{36,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(substance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}},        color={127,127,0}));
      connect(boundaryRear1.fore, substance3.rear) annotation (Line(
          points={{-58,52},{-28,52}},
          color={158,66,200},
          thickness=0.5));
      connect(substance3.fore, boundaryFore2.rear) annotation (Line(
          points={{-8,52},{30,52}},
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
       annotation (
        Icon(graphics,
             coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Documentation(info="<html>
<u>Tests for the rear and fore boundary.</u>
<u><br>Owner: <a href=\"mailto:marek@matfyz.cz\">Marek Matejak</a></u>
</html>"));
    end TestSubstance;

    model TestExternalSubstance
       extends Modelica.Icons.Example;
      Chemical.Solution solution(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL                                                        "Ideal Gas from MSL")
                                 annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

      replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
        Chemical.Interfaces.StateOfMatter "Gas substance model"
        annotation (choices(
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      ExternalSubstance
                externalSubstance1(
        useRear=true,
        useFore=false,
        redeclare package stateOfMatter = gasModel,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        FixedValue=1)     annotation (Placement(transformation(extent={{24,14},{44,34}})));

      ExternalSubstance
                externalSubstance2(
        useRear=false,
        useFore=true,
        useSolution=true,
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=10)                          annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));

      BoundaryFore boundaryFore(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      ExternalSubstance externalIdealGas(
        useRear=false,
        useFore=true,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        FixedValue=15)                                        annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
      BoundaryFore boundaryFore1(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{32,72},{52,92}})));
      BoundaryRear boundaryRear1(
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
      ExternalSubstance
                externalSubstance(
        useRear=true,
        useFore=true,       redeclare package stateOfMatter = gasModel,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3,
        FixedValue=10)    annotation (Placement(transformation(extent={{-26,42},{-6,62}})));

      BoundaryFore boundaryFore2(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{30,42},{50,62}})));
      BoundaryRear boundaryRear2(
        redeclare package stateOfMatter = gasModel,
                                 substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      ExternalSubstance
                externalSubstance4(
        useRear=true,        useFore=false,
        useSolution=true,
        redeclare package stateOfMatter = gasModel,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=1)     annotation (Placement(transformation(extent={{34,-82},{54,-62}})));

      BoundaryRear boundaryRear3(
        redeclare package stateOfMatter = gasModel,
                                 substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      ExternalSubstance
                externalSubstance3(
        useRear=true,
        useFore=true,
        useSolution=true,
        redeclare package stateOfMatter = gasModel,
        quantity=Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg,
        FixedValue=1)     annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));

      BoundaryFore boundaryFore3(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryRear boundaryRear4(
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      Substance solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{70,-82},{90,-62}})));
    equation
      connect(externalSubstance2.fore, boundaryFore.rear) annotation (Line(
          points={{-48,-16},{36,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(externalSubstance2.solution, solution.solution) annotation (Line(points={{-64,-26},{-64,-104},{60,-104},{60,-98.94}}, color={127,127,0}));
      connect(externalIdealGas.fore, boundaryFore1.rear) annotation (Line(
          points={{-52,82},{32,82}},
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
       extends Modelica.Icons.Example;
      Chemical.Solution solution(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGasMSL                                                        "Ideal Gas from MSL")
                                 annotation (Placement(transformation(extent={{-100,-100},{100,6}})));

      replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
        Chemical.Interfaces.StateOfMatter "Gas substance model"
        annotation (choices(
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      ExternalGas
                externalGas1(
        useRear=true,
        useFore=false,
        redeclare package gasModel = gasModel,
                       PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{24,14},{44,34}})));
      ExternalGas
                externalGas2(
        useRear=false,
        useFore=true,
        useSolution=true,
        redeclare package gasModel = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        PartialPressure(displayUnit="mmHg") = 1333.22387415)
                                                annotation (Placement(transformation(extent={{-68,-26},{-48,-6}})));
      BoundaryFore boundaryFore(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{36,-26},{56,-6}})));
      ExternalGas externalIdealGas(
        useRear=false,
        useFore=true,
        redeclare package gasModel = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        PartialPressure(displayUnit="mmHg") = 1999.835811225) annotation (Placement(transformation(extent={{-72,72},{-52,92}})));
      BoundaryFore boundaryFore1(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{32,72},{52,92}})));
      BoundaryRear boundaryRear1(
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-76,42},{-56,62}})));
      ExternalGas
                externalGas(
        useRear=true,
        useFore=true,       redeclare package gasModel = gasModel,                 PartialPressure(displayUnit="mmHg") = 1333.22387415)
                          annotation (Placement(transformation(extent={{-24,42},{-4,62}})));
      BoundaryFore boundaryFore2(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{30,42},{50,62}})));
      BoundaryRear boundaryRear2(
        redeclare package stateOfMatter = gasModel,
                                 substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-66,-82},{-46,-62}})));
      ExternalGas
                externalGas4(
        useRear=true,        useFore=false,
        useSolution=true,
        redeclare package gasModel = gasModel,
        PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{34,-82},{54,-62}})));
      BoundaryRear boundaryRear3(
        redeclare package stateOfMatter = gasModel,
                                 substanceData=Chemical.Substances.IdealGasesMSL.H2O(), solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-68,-54},{-48,-34}})));
      ExternalGas
                externalGas3(
        useRear=true,
        useFore=true,
        useSolution=true,
        redeclare package gasModel = gasModel,       PartialPressure(displayUnit="mmHg") = 133.322387415)
                          annotation (Placement(transformation(extent={{-18,-54},{2,-34}})));
      BoundaryFore boundaryFore3(redeclare package stateOfMatter = gasModel)
                                annotation (Placement(transformation(extent={{40,-54},{60,-34}})));
      BoundaryRear boundaryRear4(
        redeclare package stateOfMatter = gasModel,
        substanceData=Chemical.Substances.IdealGasesMSL.H2O(),
        solutionFromInput=false)
                       annotation (Placement(transformation(extent={{-74,14},{-54,34}})));
      Substance solvent(useFore=false, useSolution=true) annotation (Placement(transformation(extent={{66,-76},{86,-56}})));
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
        potentialFromInput=true,
        u0_par=110000) annotation (Placement(transformation(extent={{22,72},{42,92}})));
      inner Chemical.DropOfCommons dropOfCommons(n_flow_reg=0.01) annotation (Placement(transformation(extent={{-88,72},{-68,92}})));
      Modelica.Blocks.Sources.Step step(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{62,76},{50,88}})));
      BoundaryFore boundary_fore1(
        potentialFromInput=true,
        u0_par=110000) annotation (Placement(transformation(extent={{22,46},{42,66}})));
      Modelica.Blocks.Sources.Step step1(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{62,50},{50,62}})));
      BoundaryRear boundary_rear1(
        potentialFromInput=true,
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
        solutionFromInput=true,
        u0_par=100000,
        fore(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-22})));
      BoundaryFore boundary_fore2(potentialFromInput=true, u0_par=110000)
                       annotation (Placement(transformation(extent={{-10,-32},{10,-12}})));
      Modelica.Blocks.Sources.Step step3(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{30,-28},{18,-16}})));
      BoundaryFore boundary_fore3(potentialFromInput=true, u0_par=110000)
                       annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
      Modelica.Blocks.Sources.Step step4(
        height=-100000,
        offset=140000,
        startTime=5)
        annotation (Placement(transformation(extent={{64,-56},{52,-44}})));
      BoundaryRear boundary_rear3(
        substanceData=Chemical.Substances.Water_liquid(),
        solutionFromInput=true,
        potentialFromInput=true,
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
      TerminalInflow terminalInflow(useSolution=false) annotation (Placement(transformation(extent={{-48,50},{-28,70}})));
      TerminalOutflow terminalOutflow1(SubstanceFlow=0) annotation (Placement(transformation(extent={{66,-86},{86,-66}})));
      TerminalInflow terminalInflow1 annotation (Placement(transformation(extent={{-48,-54},{-28,-34}})));
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

    partial model PartialBoundary

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

      parameter stateOfMatter.SubstanceDataParameters substanceData    "Definition of the substance"
        annotation (choicesAllMatching = true, Dialog(enable = not useRear));

      outer Modelica.Fluid.System system "System wide properties";

      parameter Boolean useRear = false "Use rearwards conector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));
      parameter Boolean useFore = true "Use forwards connector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L
       annotation( Dialog(tab = "Advanced"));

      parameter Modelica.Units.SI.MolarFlowRate n_flow_reg=dropOfCommons.n_flow_reg "Regularization threshold of mass flow rate"
        annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Rear rear(
        redeclare package stateOfMatter = stateOfMatter,
        n_flow=n_flow_rear,
        r=r_rear_port,
        state_rearwards=state_out,
        solution=solutionState) if useRear
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

      Chemical.Interfaces.Fore fore(
        redeclare package stateOfMatter = stateOfMatter,
        n_flow=n_flow_fore,
        r=r_fore_port,
        state_forwards=state_out,
        solution=solutionState,
        definition=substanceDataVar) if useFore
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

    //  stateOfMatter.SubstanceData substanceDataVar;

      Modelica.Units.SI.MolarFlowRate n_flow "Molar change of the amount of base substance";
      Modelica.Units.SI.EnthalpyFlowRate h_flow "Change of enthalpy";

      Modelica.Units.SI.AmountOfSubstance n
        "Amount of base molecules inside all clusters in compartment";

    protected

      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

      outer Chemical.DropOfCommons dropOfCommons;

      stateOfMatter.InputSubstanceData substanceDataVar;//substanceDefinition;
      Chemical.Interfaces.SolutionState solutionState;

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

      Chemical.Interfaces.InputSubstanceState state_in_rear;
      Chemical.Interfaces.InputSubstanceState state_in_fore;
      Chemical.Interfaces.SubstanceState state_out;

    equation

      der(n_flow_rear)*L = r_rear_port - r_rear_intern;
      der(n_flow_fore)*L = r_fore_port - r_fore_intern;

     // connect(rear.definition, substanceDefinition);
      connect(rear.definition, substanceDataVar);
     // substanceDataVar = substanceDefinition;
      connect(state_in_rear,rear.state_forwards);
      connect(state_in_fore,fore.state_rearwards);

      if not useRear then
        r_rear_port = 0;
        n_flow_rear = 0;
        state_in_rear.h = 0;
        substanceDataVar = substanceData;
      end if;

      if not useFore then
        r_fore_port = 0;
        n_flow_fore = 0;
        state_in_fore.h = 0;
      end if;

      n_flow = n_flow_rear + n_flow_fore;
      h_flow = h_flow_rear + h_flow_fore;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialBoundary;

    partial model PartialSubstance
      extends PartialBoundary;

      stateOfMatter.BaseProperties substance(
       substanceDataVar=substanceDataVar,
       solutionState=solutionState,
       FixedSubstanceData=not useRear,
       substanceData=substanceData,
       amountOfBaseMolecules=n,
       m_start=m_start,
       n_flow=n_flow,
       h_flow=h_flow);

    parameter Boolean useSolution = false "Use solution connector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
        annotation (Dialog(enable=not useSolution and not useRear));

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

    protected

         Chemical.Interfaces.SolutionState solutionPortState;

    equation

      state_out.u = substance.u;
      state_out.h = substance.h;

      if (useSolution and not useRear) or (not useSolution) then
        solutionState=solutionPortState;
      end if;

      if not useSolution and not useRear then
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

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialSubstance;

    partial model PartialTerminalRear

     outer Modelica.Fluid.System system "System wide properties";

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

      Chemical.Interfaces.Rear rear(redeclare package stateOfMatter = stateOfMatter) "The substance"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

      stateOfMatter.BaseProperties substance(
        SolutionObserverOnly=true,
        substanceDataVar=rear.definition,
        solutionState=rear.solution,
        FixedSubstanceData=false,
        m_start=1,
        n_flow=0,
        h_flow=0);

      parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
      parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    protected
      Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

    initial equation
      u = u_start;
    equation

     substance.u = rear.state_forwards.u;

      TC * der(u) = rear.r;
      rear.state_rearwards.u = u;
      rear.state_rearwards.h = rear.stateOfMatter.molarEnthalpy(rear.definition,rear.solution.T);

     annotation (
       Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialTerminalRear;

    model PartialSolutionSensor

      parameter Boolean useSolution = true "Use solution connector?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
        annotation (Dialog(enable=not useSolution and not useRear));

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

    protected
        parameter Boolean useSolutionFromRear = false "Use solution from Rear port?";

        Chemical.Interfaces.SolutionState solutionPortState;

    equation

     if (useSolution and not useSolutionFromRear) or (not useSolution) then
        solutionState=solutionPortState;
      end if;

      if not useSolution and not useSolutionFromRear then
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

       annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));

    end PartialSolutionSensor;

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
