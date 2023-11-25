within Chemical;
package Boundaries "This package contains boundary models for the stream."
extends Modelica.Icons.SourcesPackage;

  model Substance "Substance in solution"
    extends Icons.Substance;

    Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
      "Molar concentration of particles";

    extends Internal.PartialSubstanceInSolutionWithAdditionalPorts;

    parameter Boolean use_mass_start = true "use mass_start, otherwise amountOfSubstance_start"
      annotation (Evaluate=true, choices(checkBox=true), Dialog(group="Initialization"));

  parameter Modelica.Units.SI.Mass mass_start=1
    "Initial mass of the substance"
    annotation (HideResult=not use_mass_start, Dialog(group="Initialization", enable=use_mass_start));

  parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=1
    "Initial amount of substance base molecules"
      annotation (HideResult=use_mass_start, Dialog(group="Initialization", enable=not use_mass_start));

    Modelica.Units.SI.Mass mass=amountOfBaseMolecules*
        molarMassOfBaseMolecule "Mass";

    parameter Boolean calculateClusteringHeat = true "Only for self clustering substances"
        annotation(Evaluate=true, choices(checkBox=true), Dialog(tab = "Clustering", enable = stateOfMatter.selfClustering(substanceData)));

  protected
    parameter Modelica.Units.SI.Mass m_start=if use_mass_start then mass_start else
      amountOfSubstance_start*molarMassOfBaseMolecule;

    parameter Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceData);

    Modelica.Units.SI.AmountOfSubstance amountOfBaseMolecules(start=
         m_start/molarMassOfBaseMolecule)
      "Amount of base molecules inside all clusters in compartment";

    Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
         m_start*stateOfMatter.specificAmountOfFreeBaseMolecule(
                                     substanceData,
                                     T=system.T_ambient,
                                     p=system.p_ambient))
      "Amount of free molecules not included inside any clusters in compartment";

    Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
         m_start*stateOfMatter.specificAmountOfParticles(
                                     substanceData,
                                     T=system.T_ambient,
                                     p=system.p_ambient))
      "Amount of particles/clusters in compartment";

    Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
        Modelica.Constants.R*solution.T))
      "Dissociation constant of hydrogen bond between base molecules";

    Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
        stateOfMatter.selfClusteringBondEnthalpy(substanceData)
      - solution.T * stateOfMatter.selfClusteringBondEntropy(substanceData)
      "Gibbs energy of hydrogen bond between H2O molecules";

    Modelica.Units.SI.AmountOfSubstance amountOfBonds
      "Amount of hydrogen bonds between molecules in compartment";

    Real logn(stateSelect=StateSelect.prefer, start=log(m_start/molarMassOfBaseMolecule))
    "Natural logarithm of the amount of base molecules in solution";

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  initial equation

    amountOfBaseMolecules = m_start/molarMassOfBaseMolecule;

  equation

    if stateOfMatter.selfClustering(substanceData) then

      //Liquid cluster theory - equilibrium:
      //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
      //amountOfParticles/solution.n = x/(1-K*x);                //sum(x[i])
      //amountOfBaseMolecules/solution.n = x/((1-K*x)^2);            //sum(i*x[i])
      //amountOfHydrogenBonds/solution.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

      amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

      //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
      x = ((2*SelfClustering_K+solution.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solution.n/amountOfBaseMolecules)+(solution.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

      amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

      //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
      //TODO: more precise calculation of other properties

     //der(enthalpy) = solution.dH + q*actualStream(port_a.h_outflow);
     //enthalpy = molarEnthalpy*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
      solution.dH =if (EnthalpyNotUsed) then 0 else der(molarEnthalpy)*
        amountOfBaseMolecules + q*molarEnthalpy - n_flow_out*h_out - n_flow_in*h_in + (
        if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
        substanceData)*der(amountOfBonds) else 0)
                      "heat transfer from other substances in solution [J/s]";

      solution.Gj =amountOfBaseMolecules*u_out + amountOfBonds*SelfClustering_dG
                      "Gibbs energy of the substance";

    else

      amountOfParticles = amountOfFreeMolecule;
      amountOfBaseMolecules = amountOfFreeMolecule;
      amountOfBonds = 0;

      //der(enthalpy) = solution.dH + q*actualStream(port_a.h_outflow);
      //enthalpy = molarEnthalpy*amountOfBaseMolecules;
      solution.dH =
        if (EnthalpyNotUsed) then  0
        else    der(molarEnthalpy)*amountOfBaseMolecules + q*molarEnthalpy
                - n_flow_out*h_out - n_flow_in*h_in
                "heat transfer from other substances in solution [J/s]";

      solution.Gj = amountOfBaseMolecules*u_out "Gibbs energy of the substance [J]";

    end if;

    //The main accumulation equation is "der(amountOfBaseMolecules)=q"
    // However, the numerical solvers can handle it in form of log(n) much better. :-)
    der(logn) = (q/amountOfBaseMolecules) "accumulation of amountOfBaseMolecules=exp(logn) [mol]";
    amountOfBaseMolecules = exp(logn);

    x = amountOfFreeMolecule/solution.n "mole fraction [mol/mol]";

    c = amountOfParticles/solution.V "concentration [mol/m3]";

    //solution flows
    solution.i = Modelica.Constants.F*z*q +
        Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
    solution.dV = molarVolume*q + der(molarVolume)*amountOfBaseMolecules "change of substance volume [m3/s]";

    //extensive properties
    solution.nj = amountOfParticles;
    solution.mj = amountOfBaseMolecules*molarMassOfBaseMolecule;
    solution.Vj = amountOfBaseMolecules*molarVolume;
    solution.Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
    solution.Ij = (1/2)*(amountOfBaseMolecules*z^2);

       annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={128,0,255},
            textString="%name")}), Documentation(revisions="<html>
<p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
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

  model Source "Boundary model of a source"


    parameter Boolean potentialFromInput = false "Use input connector for electro-chemical potential?";
    parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";

    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "Electro-chemical potential set value" annotation (Dialog(enable=not potentialFromInput));
    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "Molar enthalpy set value"
      annotation (Dialog(enable=not enthalpyFromInput));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J.mol-1") if potentialFromInput "Electro-chemical potential input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,40},{0,80}}), iconTransformation(extent={{-40,40},{0,80}})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/mol") if enthalpyFromInput "Enthalpy input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,-40},{0,0}}), iconTransformation(extent={{-40,-20},{0,20}})));
    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(extent={{80,-20},{120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal electro-chemical potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

  equation

     connect(u0_var, u0);
     if not potentialFromInput then
       u0 = u0_par;
     end if;

     connect(h0_var, h0);
     if not enthalpyFromInput then
       h0 = h0_par;
     end if;

    L*der(outlet.n_flow) = outlet.r - 0;
    outlet.u = u0;
    outlet.h = h0;

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
            points={{60,0},{100,0}},
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
    <p>Source of a Chemical substance stream. The state can be given as fix values or as a real signal. </p>
<p>Before its inertance the source has an inertial potential of 0 by definition.</p>
</html>"));
  end Source;

  model Sink "Boundary model of sink"

    parameter Boolean potentialFromInput = false "If true electro-chemical potential comes from real input";
    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "Electro-chemical potential setpoint of Sink" annotation (Dialog(enable=not pressureFromInput));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of electro-chemical potential" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol") if potentialFromInput "Potential setpoint [J/mol]"
      annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={20,0}), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={20,0})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal electro-chemical potential connector";
    Modelica.Units.SI.ChemicalPotential r;

    Modelica.Units.SI.ChemicalPotential u=inlet.u;

  equation

    connect(u0_var, u0);
    if not potentialFromInput then
      u0 = u0_par;
    end if;

    der(inlet.n_flow)*L = inlet.r - r;
    r + u = u0;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-56,76},{4,-84}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-60,80},{0,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-100,0},{-60,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{-60,80},{-60,-80}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{-12,80},{-12,-80}},
            color={255,255,255},
            thickness=1),
          Line(
            points={{-28,80},{-28,-80}},
            color={255,255,255},
            thickness=0.5),
          Line(points={{-44,80},{-44,-80}}, color={255,255,255})}), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Sink for a thermofluid stream. The pressure can be set or given by a real signal via input connector.</p>
<p>The inertial pressure after the sinks inertance is by definition the difference between the input pressure and the set pressure. The sink therefore acts by definition as the origin of the energy to accelerate the stream. </p>
</html>"));
  end Sink;

  model TerminalSource "Source that imposes n_flow = 0"

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarEnthalpy h=0 "Source enthalpy";
    parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial potential";

    Chemical.Interfaces.Outlet outlet
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_0;

  equation
    outlet.n_flow = 0;

    TC * der(u) = outlet.r;
    outlet.u = u;
    outlet.h = h;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{34,26},{74,-34}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{50,0},{100,0}},
            color={158,66,200},
            thickness=0.5),
          Rectangle(
            extent={{30,30},{70,-30}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{70,30},{30,-30}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{30,30},{70,-30}},
            color={158,66,200},
            thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Source that terminates the flow. </p>
<p>It imposes a n_flow=0 boundary and with a time constant, adapts the pressure such that the inertial electro-chemical potetntial r goes to zero.</p>
</html>"));
  end TerminalSource;

  model TerminalSink "Sink that imposes m_flow=0"

    Chemical.Interfaces.Inlet inlet
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

  equation
    inlet.n_flow = 0;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Line(
            points={{-100,0},{-50,0}},
            color={158,66,200},
            thickness=0.5),
          Rectangle(
            extent={{-56,26},{-16,-34}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-60,30},{-20,-30}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{-20,30},{-60,-30}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{-60,30},{-20,-30}},
            color={158,66,200},
            thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Sink that terminates the flow. </p>
<p>It imposes a m_flow=0 boundary.</p>
</html>"));
  end TerminalSink;

  package Tests "Tests for the boundary package"
    extends Modelica.Icons.ExamplesPackage;

    model SourceSink "Test for source and sink model"
      extends Modelica.Icons.Example;

      inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{58,-74},{78,-54}})));
      Source source(
        u0_par=-200000, L=0,
        outlet(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-32,10},{-12,30}})));
      Sink sink(u0_par=-100000)
        annotation (Placement(transformation(extent={{14,10},{34,30}})));
      Source source1(
        u0_par=-200000,
        outlet(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-32,-30},{-12,-10}})));
      Sink sink1(
        L=0,
        u0_par=-100000)
        annotation (Placement(transformation(extent={{14,-30},{34,-10}})));
    equation
      connect(sink.inlet, source.outlet) annotation (Line(
          points={{14,20},{-12,20}},
          color={28,108,200},
          thickness=0.5));
      connect(sink1.inlet, source1.outlet) annotation (Line(
          points={{14,-20},{-12,-20}},
          color={28,108,200},
          thickness=0.5));

      annotation (experiment(StopTime=1, Tolerance=1e-6, Interval=0.001),
        Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end SourceSink;

    model Volumes "Test Volumes"
      extends Modelica.Icons.Example;

      replaceable package Medium = ThermofluidStream.Media.myMedia.Air.SimpleAir constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
        "Medium package" annotation (Documentation(info="<html>
<p>
Medium package used in the Test.
</p>
</html>"));

      package MediumMix = ThermofluidStream.Media.myMedia.IdealGases.MixtureGases.CombustionAir
        "Medium package"
        annotation (Documentation(info="<html>
<p>
Medium package used in the Test of the MixVolumes.
</p>
</html>"));

      inner ThermofluidStream.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
        annotation (Placement(transformation(extent={{70,-170},{90,-150}})));
      Source source(redeclare package Medium = Medium,
        p0_par=200000,
        outlet(m_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
      Sink sink(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{50,-60},{70,-40}})));
      Source source1(redeclare package Medium = Medium, p0_par=200000)
        annotation (Placement(transformation(extent={{-70,-90},{-50,-70}})));
      Sink sink1(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{50,-90},{70,-70}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
            displayUnit="K") = 500)
        annotation (Placement(transformation(extent={{-48,-120},{-28,-100}})));
      ThermofluidStream.Processes.FlowResistance flowResistance(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{20,-90},{40,-70}})));
      ThermofluidStream.Processes.FlowResistance flowResistance1(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
      Source source2(redeclare package Medium = MediumMix, p0_par=200000)
        annotation (Placement(transformation(extent={{-70,10},{-50,30}})));
      Source source3(redeclare package Medium = MediumMix, p0_par=200000)
        annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
      Sink sink2(redeclare package Medium = MediumMix, p0_par=100000)
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));
      ThermofluidStream.Processes.FlowResistance flowResistance2(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
      ThermofluidStream.Processes.FlowResistance flowResistance3(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      ThermofluidStream.Processes.FlowResistance flowResistance4(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.01,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      Source source4(
        redeclare package Medium = MediumMix,
        p0_par=200000,
        Xi0_par={1,0})
        annotation (Placement(transformation(extent={{-70,80},{-50,100}})));
      Source source5(
        redeclare package Medium = MediumMix,
        p0_par=200000,
        Xi0_par={0,1})
        annotation (Placement(transformation(extent={{-70,40},{-50,60}})));
      Sink sink3(redeclare package Medium = MediumMix, p0_par=100000)
        annotation (Placement(transformation(extent={{50,60},{70,80}})));
      ThermofluidStream.Processes.FlowResistance flowResistance5(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e5))
        annotation (Placement(transformation(extent={{-40,80},{-20,100}})));
      ThermofluidStream.Processes.FlowResistance flowResistance6(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=2e5))
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      ThermofluidStream.Processes.FlowResistance flowResistance8(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e1))
        annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Source source6(redeclare package Medium = Medium, p0_par=180000)
        annotation (Placement(transformation(extent={{-70,-160},{-50,-140}})));
      Sink sink4(redeclare package Medium = Medium, p0_par=130000)
        annotation (Placement(transformation(extent={{50,-200},{70,-180}})));
      ThermofluidStream.Processes.FlowResistance flowResistance7(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=1,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{-40,-160},{-20,-140}})));
      ThermofluidStream.Processes.FlowResistance flowResistance9(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarTurbulentPressureLoss (material=ThermofluidStream.Processes.Internal.Material.wood))
        annotation (Placement(transformation(extent={{20,-200},{40,-180}})));
      Source source7(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{-70,110},{-50,130}})));
      Sink sink5(redeclare package Medium = Medium, p0_par=200000)
        annotation (Placement(transformation(extent={{50,110},{70,130}})));
      ThermofluidStream.Processes.FlowResistance flowResistance10(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="mm") = 0.01,
        l=1,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e5))
        annotation (Placement(transformation(extent={{20,110},{40,130}})));
      ThermofluidStream.Processes.FlowResistance flowResistance11(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="mm") = 0.01,
        l=1,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e5))
        annotation (Placement(transformation(extent={{-40,110},{-20,130}})));
      Source source8(
        redeclare package Medium = MediumMix,
        p0_par=200000,
        Xi0_par={1,0})
        annotation (Placement(transformation(extent={{-70,180},{-50,200}})));
      Source source9(
        redeclare package Medium = MediumMix,
        p0_par=100000,
        Xi0_par={0,1})
        annotation (Placement(transformation(extent={{-70,140},{-50,160}})));
      Sink sink6(redeclare package Medium = MediumMix, p0_par=300000)
        annotation (Placement(transformation(extent={{50,160},{70,180}})));
      ThermofluidStream.Processes.FlowResistance flowResistance12(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e5))
        annotation (Placement(transformation(extent={{-40,180},{-20,200}})));
      ThermofluidStream.Processes.FlowResistance flowResistance13(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=2e5))
        annotation (Placement(transformation(extent={{-40,140},{-20,160}})));
      ThermofluidStream.Processes.FlowResistance flowResistance14(
        redeclare package Medium = MediumMix,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r=0.1,
        l=10,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=1e1))
        annotation (Placement(transformation(extent={{20,160},{40,180}})));
    equation
      connect(volume.inlet, source.outlet) annotation (Line(
          points={{-30,-50},{-50,-50}},
          color={28,108,200},
          thickness=0.5));
      connect(volume.outlet, flexVolume1.inlet) annotation (Line(
          points={{-10,-50},{10,-50}},
          color={28,108,200},
          thickness=0.5));
      connect(flexVolume1.outlet, sink.inlet) annotation (Line(
          points={{30,-50},{50,-50}},
          color={28,108,200},
          thickness=0.5));

      connect(fixedTemperature.port,heatportVolume. heatPort)
        annotation (Line(points={{-28,-110},{0,-110},{0,-88}}, color={191,0,0}));
      connect(heatportVolume.outlet,flowResistance. inlet) annotation (Line(
          points={{10,-80},{20,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance.outlet,sink1. inlet) annotation (Line(
          points={{40,-80},{50,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance1.inlet,source1. outlet) annotation (Line(
          points={{-40,-80},{-50,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance1.outlet,heatportVolume. inlet) annotation (Line(
          points={{-20,-80},{-10,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(source2.outlet, flowResistance2.inlet) annotation (Line(
          points={{-50,20},{-40,20}},
          color={28,108,200},
          thickness=0.5));
      connect(source3.outlet, flowResistance3.inlet) annotation (Line(
          points={{-50,-20},{-40,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(volumeMix.outlet, flowResistance4.inlet) annotation (Line(
          points={{10,0},{20,0}},
          color={28,108,200},
          thickness=0.5));
      connect(sink2.inlet, flowResistance4.outlet) annotation (Line(
          points={{50,0},{40,0}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance3.outlet, volumeMix.inlet[1]) annotation (Line(
          points={{-20,-20},{-16,-20},{-16,-0.5},{-10,-0.5}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance2.outlet, volumeMix.inlet[2]) annotation (Line(
          points={{-20,20},{-16,20},{-16,0.5},{-10,0.5}},
          color={28,108,200},
          thickness=0.5));
      connect(source4.outlet,flowResistance5. inlet) annotation (Line(
          points={{-50,90},{-40,90}},
          color={28,108,200},
          thickness=0.5));
      connect(source5.outlet,flowResistance6. inlet) annotation (Line(
          points={{-50,50},{-40,50}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance6.outlet, volumeMix1.inlet[1]) annotation (Line(
          points={{-20,50},{-16,50},{-16,69.5},{-10,69.5}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance5.outlet, volumeMix1.inlet[2]) annotation (Line(
          points={{-20,90},{-16,90},{-16,70.5},{-10,70.5}},
          color={28,108,200},
          thickness=0.5));
      connect(sink3.inlet, flowResistance8.outlet) annotation (Line(
          points={{50,70},{40,70}},
          color={28,108,200},
          thickness=0.5));
      connect(volumeMix1.outlet, flowResistance8.inlet) annotation (Line(
          points={{10,70},{20,70}},
          color={28,108,200},
          thickness=0.5));
      connect(source6.outlet, flowResistance7.inlet) annotation (Line(
          points={{-50,-150},{-40,-150}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume1.inlet, flowResistance7.outlet)
        annotation (Line(
          points={{-8,-150},{-20,-150}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume2.outlet, flowResistance9.inlet)
        annotation (Line(
          points={{10,-190},{20,-190}},
          color={28,108,200},
          thickness=0.5));
      connect(sink4.inlet, flowResistance9.outlet) annotation (Line(
          points={{50,-190},{40,-190}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume3.heatPort, heatportVolume.heatPort) annotation (Line(points={{42,-110},{0,-110},{0,-88}},
                                                                                                                     color={191,0,0}));
      connect(volumeMix2.heatPort, heatportVolume.heatPort) annotation (Line(points={{42,-140},{28,-140},{28,-110},{0,-110},{0,-88}},
                                                                                                                                  color={191,0,0}));
      connect(volume1.outlet, flowResistance10.inlet) annotation (Line(
          points={{8,120},{20,120}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance10.outlet, sink5.inlet) annotation (Line(
          points={{40,120},{50,120}},
          color={28,108,200},
          thickness=0.5));
      connect(volume1.inlet, flowResistance11.outlet) annotation (Line(
          points={{-12,120},{-20,120}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance11.inlet, source7.outlet) annotation (Line(
          points={{-40,120},{-50,120}},
          color={28,108,200},
          thickness=0.5));
      connect(source8.outlet, flowResistance12.inlet) annotation (Line(
          points={{-50,190},{-40,190}},
          color={28,108,200},
          thickness=0.5));
      connect(source9.outlet, flowResistance13.inlet) annotation (Line(
          points={{-50,150},{-40,150}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance13.outlet, volumeMix3.inlet[1])
        annotation (Line(
          points={{-20,150},{-14,150},{-14,169.5},{-8,169.5}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance12.outlet, volumeMix3.inlet[2])
        annotation (Line(
          points={{-20,190},{-14,190},{-14,170.5},{-8,170.5}},
          color={28,108,200},
          thickness=0.5));
      connect(sink6.inlet, flowResistance14.outlet) annotation (Line(
          points={{50,170},{40,170}},
          color={28,108,200},
          thickness=0.5));
      connect(volumeMix3.outlet, flowResistance14.inlet) annotation (Line(
          points={{12,170},{20,170}},
          color={28,108,200},
          thickness=0.5));
      annotation (experiment(StopTime=1, Tolerance=1e-6, Interval=0.001),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-200},{100,200}})),
            Documentation(info="<html>
<p>This test purposely violates the asserts for reversed mass-flows to test the behavior of volumes for prolonged reversed massflow for some of the volumes. The assert is set to warning.</p>
<p><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end Volumes;

    model TerminalSourceSink "Test for Terminal source and sink model"
      extends Modelica.Icons.Example;

      replaceable package Medium = ThermofluidStream.Media.myMedia.Air.SimpleAir constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
        "Medium package" annotation (Documentation(info="<html>
<p>
Medium package used in the Test.
</p>
</html>"));

      inner ThermofluidStream.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{50,-70},{70,-50}})));
      Sink sink2(redeclare package Medium = Medium,
        pressureFromInput=true,
        p0_par=100000)
        annotation (Placement(transformation(extent={{8,4},{28,24}})));
      Modelica.Blocks.Sources.Pulse pulse(
        amplitude=5e4,
        width=35,
        period=2,
        offset=0.8e5,
        startTime=0.2)
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      TerminalSource terminalSource(redeclare package Medium = Medium, TC=0.1) annotation (Placement(transformation(extent={{-26,4},{-6,24}})));
      TerminalSink terminalSink(redeclare package Medium = Medium) annotation (Placement(transformation(extent={{8,-26},{28,-6}})));
      Source source(redeclare package Medium = Medium, pressureFromInput=true) annotation (Placement(transformation(extent={{-26,-26},{-6,-6}})));
    equation

      connect(pulse.y, sink2.p0_var) annotation (Line(points={{-39,0},{26,0},{26,14},{20,14}}, color={0,0,127}));
      connect(sink2.inlet, terminalSource.outlet) annotation (Line(
          points={{8,14},{-6,14}},
          color={28,108,200},
          thickness=0.5));
      connect(source.outlet, terminalSink.inlet) annotation (Line(
          points={{-6,-16},{8,-16}},
          color={28,108,200},
          thickness=0.5));
      connect(source.p0_var, sink2.p0_var) annotation (Line(points={{-18,-10},{-26,-10},{-26,0},{26,0},{26,14},{20,14}}, color={0,0,127}));
      annotation (experiment(StopTime=10, Tolerance=1e-6, Interval=0.01),
        Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end TerminalSourceSink;

    model VolumesDirectCoupling "Test Volumes"
      extends Modelica.Icons.Example;

      replaceable package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater constrainedby ThermofluidStream.Media.myMedia.Interfaces.PartialMedium
        "Medium package"
        annotation (choicesAllMatching=true, Documentation(info="<html>
<p>
Medium package used in the Test.
</p>
</html>"));

      package MediumMix = ThermofluidStream.Media.myMedia.IdealGases.MixtureGases.CombustionAir
        "Medium package"
        annotation (Documentation(info="<html>
<p>
Medium package used in the Test of the MixVolumes.
</p>
</html>"));

      inner ThermofluidStream.DropOfCommons dropOfCommons(k_volume_damping=0.5, assertionLevel=AssertionLevel.warning)
        annotation (Placement(transformation(extent={{70,-148},{90,-128}})));
      Sink sink4(redeclare package Medium = Medium, p0_par=130000)
        annotation (Placement(transformation(extent={{36,-10},{56,10}})));
      Source source(redeclare package Medium = Medium,
        p0_par=200000,
        outlet(m_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-58,-96},{-38,-76}})));
      Sink sink1(redeclare package Medium = Medium, p0_par=130000)
        annotation (Placement(transformation(extent={{36,-30},{56,-10}})));
      Source source1(redeclare package Medium = Medium,
        p0_par=200000,
        outlet(m_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-58,-118},{-38,-98}})));
      Sink sink2(redeclare package Medium = Medium, p0_par=130000)
        annotation (Placement(transformation(extent={{36,-50},{56,-30}})));
      Sink sink3(redeclare package Medium = Medium, p0_par=101000)
        annotation (Placement(transformation(extent={{36,-76},{56,-56}})));
    equation

      connect(heatportVolume2.outlet, sink4.inlet) annotation (Line(
          points={{10,0},{36,0}},
          color={28,108,200},
          thickness=0.5));
      connect(source.outlet, heatportVolume1.inlet) annotation (Line(
          points={{-38,-86},{-10,-86}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume4.inlet, heatportVolume3.outlet) annotation (Line(
          points={{10,30},{-10,30}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume6.inlet, heatportVolume5.outlet) annotation (Line(
          points={{10,60},{-10,60}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume7.outlet, sink1.inlet) annotation (Line(
          points={{10,-20},{36,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume9.inlet, heatportVolume8.outlet) annotation (Line(
          points={{10,88},{-10,88}},
          color={28,108,200},
          thickness=0.5));
      connect(source1.outlet, heatportVolume10.inlet) annotation (Line(
          points={{-38,-108},{-10,-108}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume11.outlet, sink2.inlet) annotation (Line(
          points={{10,-40},{36,-40}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume13.inlet, heatportVolume12.outlet)
        annotation (Line(
          points={{10,112},{-10,112}},
          color={28,108,200},
          thickness=0.5));
      connect(heatportVolume14.outlet, sink3.inlet) annotation (Line(
          points={{10,-66},{36,-66}},
          color={28,108,200},
          thickness=0.5));
      annotation (experiment(
          StopTime=0.05,
          Interval=5e-06,
          Tolerance=1e-07,
          __Dymola_Algorithm="Dassl"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-160},{100,160}})),
            Documentation(info="<html>
<p>This test uses deceased tolerance and saves more points (smaller output interval) because of the different time-scales and relative stiff equation system.</p>
<p><br>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end VolumesDirectCoupling;

    model DynamicBoundaries "Test for DynamicInflow and Outflow"
      extends Modelica.Icons.Example;

      package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater;

      Source source(redeclare package Medium=Medium, p0_par=101000)
        annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      DynamicPressureInflow dynamicPressureInflow(redeclare package Medium=Medium,
        areaFromInput=true,
        velocityFromInput=true)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Sink sink(redeclare package Medium=Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,-10},{80,10}})));
      ThermofluidStream.Processes.FlowResistance flowResistance(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="mm") = 0.005,
        l=10,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Source source1(redeclare package Medium=Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{-80,-50},{-60,-30}})));
      DynamicPressureInflow dynamicPressureInflow1(redeclare package Medium=Medium,
        A_par=0.005,
        v_in_par=-10,
        assumeConstantDensity=true)
        annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
      Sink sink1(redeclare package Medium=Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,-50},{80,-30}})));
      ThermofluidStream.Processes.FlowResistance flowResistance1(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        m_flow_0=-1,
        r(displayUnit="mm") = 0.005,
        l=10,
        L_value=1000,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      Modelica.Blocks.Sources.Ramp ramp(
        height=1,
        duration=0.5,
        startTime=0.4) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      inner ThermofluidStream.DropOfCommons dropOfCommons(m_flow_reg=0.01) annotation (Placement(transformation(extent={{60,20},{80,40}})));
      Modelica.Blocks.Sources.Ramp ramp1(
        height=-0.98e-3,
        duration=0.4,
        offset=1e-3,
        startTime=0) annotation (Placement(transformation(extent={{0,20},{-20,40}})));
      Source source2(redeclare package Medium = Medium, p0_par=101000)
        annotation (Placement(transformation(extent={{-80,50},{-60,70}})));
      DynamicPressureInflow dynamicPressureInflow2(
        redeclare package Medium = Medium,
        areaFromInput=false,
        velocityFromInput=false,
        assumeConstantDensity=true)
        annotation (Placement(transformation(extent={{-50,50},{-30,70}})));
      Sink sink2(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,50},{80,70}})));
      ThermofluidStream.Processes.FlowResistance flowResistance2(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="mm") = 0.005,
        l=10,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,50},{10,70}})));
      Modelica.Blocks.Sources.Ramp ramp2(
        height=1,
        duration=0.5,
        startTime=0.4) annotation (Placement(transformation(extent={{0,80},{20,100}})));
      Modelica.Blocks.Sources.Ramp ramp3(
        height=-0.95e-3,
        duration=0.4,
        offset=1e-3,
        startTime=0) annotation (Placement(transformation(extent={{80,80},{60,100}})));
      Source source3(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{-80,-70},{-60,-50}})));
      DynamicPressureInflow dynamicPressureInflow3(
        redeclare package Medium = Medium,
        A_par=0.01,
        v_in_par=0,
        assumeConstantDensity=true)
        annotation (Placement(transformation(extent={{-50,-70},{-30,-50}})));
      Sink sink3(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,-70},{80,-50}})));
      ThermofluidStream.Processes.FlowResistance flowResistance3(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        m_flow_0=-1,
        r(displayUnit="mm") = 0.005,
        l=10,
        L_value=1000,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
      Source source4(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
      DynamicPressureInflow dynamicPressureInflow4(
        redeclare package Medium = Medium,
        A_par=0.01,
        v_in_par=0,
        assumeConstantDensity=true)
        annotation (Placement(transformation(extent={{-50,-90},{-30,-70}})));
      Sink sink4(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,-90},{80,-70}})));
      ThermofluidStream.Processes.FlowResistance flowResistance4(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        m_flow_0=1,
        r(displayUnit="mm") = 0.005,
        l=10,
        L_value=1000,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,-90},{10,-70}})));
      Source source5(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
      DynamicPressureInflow dynamicPressureInflow5(
        redeclare package Medium = Medium,
        A_par=0.005,
        v_in_par=10,
        assumeConstantDensity=true)
        annotation (Placement(transformation(extent={{-50,-30},{-30,-10}})));
      Sink sink5(redeclare package Medium = Medium, p0_par=100000)
        annotation (Placement(transformation(extent={{60,-30},{80,-10}})));
      ThermofluidStream.Processes.FlowResistance flowResistance5(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        m_flow_0=-1,
        r(displayUnit="mm") = 0.005,
        l=10,
        L_value=1000,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.laminarPressureLoss)
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
    equation
      connect(source.outlet, dynamicPressureInflow.inlet) annotation (Line(
          points={{-60,0},{-50,0}},
          color={28,108,200},
          thickness=0.5));
      connect(sink.inlet, dynamicPressureOutflow.outlet) annotation (Line(
          points={{60,0},{50,0}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow.outlet, flowResistance.inlet)
        annotation (Line(
          points={{-30,0},{-10,0}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow.inlet, flowResistance.outlet)
        annotation (Line(
          points={{30,0},{10,0}},
          color={28,108,200},
          thickness=0.5));
      connect(source1.outlet, dynamicPressureInflow1.inlet)
        annotation (Line(
          points={{-60,-40},{-50,-40}},
          color={28,108,200},
          thickness=0.5));
      connect(sink1.inlet, dynamicPressureOutflow1.outlet) annotation (Line(
          points={{60,-40},{50,-40}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow1.outlet, flowResistance1.inlet)
        annotation (Line(
          points={{-30,-40},{-10,-40}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow1.inlet, flowResistance1.outlet)
        annotation (Line(
          points={{30,-40},{10,-40}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow.v_in_var, ramp.y) annotation (Line(points={{-40,10},{-40,30},{-59,30}}, color={0,0,127}));
      connect(dynamicPressureInflow.A_var, ramp1.y) annotation (Line(points={{-34,10},{-34,30},{-21,30}}, color={0,0,127}));
      connect(source2.outlet,dynamicPressureInflow2. inlet)
        annotation (Line(
          points={{-60,60},{-50,60}},
          color={28,108,200},
          thickness=0.5));
      connect(sink2.inlet,dynamicPressureOutflow2. outlet) annotation (Line(
          points={{60,60},{50,60}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow2.outlet,flowResistance2. inlet)
        annotation (Line(
          points={{-30,60},{-10,60}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow2.inlet,flowResistance2. outlet)
        annotation (Line(
          points={{30,60},{10,60}},
          color={28,108,200},
          thickness=0.5));
      connect(ramp2.y,dynamicPressureOutflow2. v_out_var) annotation (Line(points={{21,90},{34,90},{34,70}}, color={0,0,127}));
      connect(ramp3.y,dynamicPressureOutflow2. A_var) annotation (Line(points={{59,90},{40,90},{40,70}}, color={0,0,127}));
      connect(source3.outlet,dynamicPressureInflow3. inlet)
        annotation (Line(
          points={{-60,-60},{-50,-60}},
          color={28,108,200},
          thickness=0.5));
      connect(sink3.inlet,dynamicPressureOutflow3. outlet) annotation (Line(
          points={{60,-60},{50,-60}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow3.outlet,flowResistance3. inlet)
        annotation (Line(
          points={{-30,-60},{-10,-60}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow3.inlet,flowResistance3. outlet)
        annotation (Line(
          points={{30,-60},{10,-60}},
          color={28,108,200},
          thickness=0.5));
      connect(source4.outlet,dynamicPressureInflow4. inlet)
        annotation (Line(
          points={{-60,-80},{-50,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(sink4.inlet,dynamicPressureOutflow4. outlet) annotation (Line(
          points={{60,-80},{50,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow4.outlet,flowResistance4. inlet)
        annotation (Line(
          points={{-30,-80},{-10,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow4.inlet,flowResistance4. outlet)
        annotation (Line(
          points={{30,-80},{10,-80}},
          color={28,108,200},
          thickness=0.5));
      connect(source5.outlet,dynamicPressureInflow5. inlet)
        annotation (Line(
          points={{-60,-20},{-50,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(sink5.inlet,dynamicPressureOutflow5. outlet) annotation (Line(
          points={{60,-20},{50,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureInflow5.outlet,flowResistance5. inlet)
        annotation (Line(
          points={{-30,-20},{-10,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(dynamicPressureOutflow5.inlet,flowResistance5. outlet)
        annotation (Line(
          points={{30,-20},{10,-20}},
          color={28,108,200},
          thickness=0.5));
      annotation (
        experiment(StopTime=1, Tolerance=1e-6, Interval=0.001),
        Documentation(info="<html>
<p>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end DynamicBoundaries;

    model PhaseSeperator
      extends Modelica.Icons.Example;

      package Medium = ThermofluidStream.Media.myMedia.Water.StandardWater;

      Boundaries.Source source(
        redeclare package Medium = Medium,
        setEnthalpy=true,
        enthalpyFromInput=true,
        p0_par=120000,
        h0_par=2000) annotation (Placement(transformation(extent={{-90,10},{-70,30}})));
      Boundaries.Sink sink(redeclare package Medium = Medium, p0_par=100000) annotation (Placement(transformation(extent={{76,-30},{96,-10}})));

      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-56,-38},{-36,-18}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect1(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect2(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{20,30},{40,50}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect3(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{64,-52},{84,-32}})));
      Modelica.Blocks.Sources.TimeTable
                                   timeTable(
        table=[0.0,1500e3; 24.9,1500e3; 25.1,3500e3; 49.9,3500e3; 50.1,1500e3; 74.9,1500e3; 75.1,410e3; 99.9,410e3; 100.1,1500e3; 124.9,1500e3; 1e10,
            1500e3],
        offset=0,
        startTime=0) annotation (Placement(transformation(extent={{-120,-16},{-100,4}})));
      ThermofluidStream.Processes.FlowResistance flowResistance1(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="cm") = 0.05,
        l=1,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=5000))
        annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
      ThermofluidStream.Processes.FlowResistance flowResistance2(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="cm") = 0.05,
        l=1,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=5000))
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      ThermofluidStream.Processes.FlowResistance flowResistance(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="cm") = 0.05,
        l=1,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=5000))
        annotation (Placement(transformation(extent={{30,10},{50,30}})));
      ThermofluidStream.Processes.FlowResistance flowResistance3(
        redeclare package Medium = Medium,
        initM_flow=ThermofluidStream.Utilities.Types.InitializationMethods.state,
        r(displayUnit="cm") = 0.05,
        l=1,
        computeL=false,
        redeclare function pLoss = ThermofluidStream.Processes.Internal.FlowResistance.linearQuadraticPressureLoss (k=5000))
        annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect4(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-8,-50},{12,-30}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect5(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-8,30},{12,50}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{-8,38},{12,58}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect1(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{20,38},{40,58}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect2(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{-8,-58},{12,-38}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect3(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{20,-58},{40,-38}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect4(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{64,-60},{84,-40}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect5(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{-56,-46},{-36,-26}})));
      Boundaries.Source source1(
        redeclare package Medium = Medium,
        setEnthalpy=true,
        enthalpyFromInput=true,
        p0_par=120000,
        h0_par=2000) annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect6(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{-56,20},{-36,40}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect6(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{-56,28},{-36,48}})));
      Boundaries.Sink sink1(redeclare package Medium = Medium, p0_par=100000) annotation (Placement(transformation(extent={{76,10},{96,30}})));
      ThermofluidStream.Sensors.TwoPhaseSensorSelect twoPhaseSensorSelect7(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.TwoPhaseQuantities.x_kgpkg) annotation (Placement(transformation(extent={{64,32},{84,52}})));
      ThermofluidStream.Sensors.SingleSensorSelect singleSensorSelect7(
        redeclare package Medium = Medium,
        digits=2,
        quantity=ThermofluidStream.Sensors.Internal.Types.Quantities.p_bar) annotation (Placement(transformation(extent={{64,40},{84,60}})));
      inner ThermofluidStream.DropOfCommons dropOfCommons(assertionLevel=AssertionLevel.warning)
        annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
    equation
      connect(source.h0_var, timeTable.y) annotation (Line(points={{-82,20},{-90,20},{-90,-6},{-99,-6}}, color={0,0,127}));
      connect(twoPhaseSensorSelect1.inlet, accumulator.outlet)
        annotation (Line(
          points={{20,-40},{16,-40},{16,-20},{10,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect2.inlet,receiver. outlet)
        annotation (Line(
          points={{20,40},{16,40},{16,20},{10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance1.outlet,receiver. inlet) annotation (Line(
          points={{-20,20},{-10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(accumulator.inlet, flowResistance2.outlet) annotation (Line(
          points={{-10,-20},{-20,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance3.inlet, accumulator.outlet) annotation (Line(
          points={{30,-20},{10,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(flowResistance.inlet,receiver. outlet) annotation (Line(
          points={{30,20},{10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect4.inlet, flowResistance2.outlet)
        annotation (Line(
          points={{-8,-40},{-14,-40},{-14,-20},{-20,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect5.inlet,receiver. inlet)
        annotation (Line(
          points={{-8,40},{-14,40},{-14,20},{-10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect.inlet,receiver. inlet)
        annotation (Line(
          points={{-8,48},{-14,48},{-14,20},{-10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect1.inlet,receiver. outlet)
        annotation (Line(
          points={{20,48},{16,48},{16,20},{10,20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect3.inlet, accumulator.outlet)
        annotation (Line(
          points={{20,-48},{16,-48},{16,-20},{10,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect2.inlet, flowResistance2.outlet)
        annotation (Line(
          points={{-8,-48},{-14,-48},{-14,-20},{-20,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(source.outlet, flowResistance1.inlet) annotation (Line(
          points={{-70,20},{-40,20}},
          color={28,108,200},
          thickness=0.5));
      connect(source1.outlet, flowResistance2.inlet)
        annotation (Line(
          points={{-70,-20},{-56,-20},{-56,-20},{-40,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(source1.h0_var, timeTable.y) annotation (Line(points={{-82,-20},{-90,-20},{-90,-6},{-99,-6}}, color={0,0,127}));
      connect(twoPhaseSensorSelect.inlet, flowResistance2.inlet)
        annotation (Line(
          points={{-56,-28},{-60,-28},{-60,-20},{-40,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect5.inlet, flowResistance2.inlet)
        annotation (Line(
          points={{-56,-36},{-60,-36},{-60,-20},{-40,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect6.inlet, flowResistance1.inlet)
        annotation (Line(
          points={{-56,38},{-60,38},{-60,20},{-40,20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect6.inlet, flowResistance1.inlet)
        annotation (Line(
          points={{-56,30},{-60,30},{-60,20},{-40,20}},
          color={28,108,200},
          thickness=0.5));
      connect(sink.inlet, flowResistance3.outlet) annotation (Line(
          points={{76,-20},{50,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(sink1.inlet, flowResistance.outlet) annotation (Line(
          points={{76,20},{50,20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect4.inlet, flowResistance3.outlet)
        annotation (Line(
          points={{64,-50},{60,-50},{60,-20},{50,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect3.inlet, flowResistance3.outlet)
        annotation (Line(
          points={{64,-42},{60,-42},{60,-20},{50,-20}},
          color={28,108,200},
          thickness=0.5));
      connect(twoPhaseSensorSelect7.inlet, flowResistance.outlet)
        annotation (Line(
          points={{64,42},{60,42},{60,20},{50,20}},
          color={28,108,200},
          thickness=0.5));
      connect(singleSensorSelect7.inlet, flowResistance.outlet)
        annotation (Line(
          points={{64,50},{60,50},{60,20},{50,20}},
          color={28,108,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
                                                                     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,-100},{120,100}})),
        experiment(StopTime=125, Tolerance=1e-6, Interval=0.125, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>Owner: <a href=\"mailto:michael.meissner@dlr.de\">Michael Mei&szlig;ner</a></p>
</html>"));
    end PhaseSeperator;
    annotation (Documentation(info="<html>
<p>
Test package for the Boundaries package of ThermofluidStream.
</p>
</html>"));
  end Tests;

  package Internal "Partials and Internal functions"
  extends Modelica.Icons.InternalPackage;

    partial model PartialSubstance
      import Chemical;

     outer Modelica.Fluid.System system "System wide properties";

     parameter Real L=1e-5;

      parameter Boolean useInlet = false "If true inlet is added";
      parameter Boolean useOutlet = true "If true outlet is added";

      Interfaces.Inlet inlet(
        r=r_in,
        n_flow=n_flow_in,
        u=u_in,
        h=h_in) if useInlet "The substance entering"
        annotation (Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

      Chemical.Interfaces.Outlet outlet(r=r_out,n_flow=n_flow_out,u=u_out,h=h_out) if useOutlet "The substance exiting"
                                                        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

     replaceable package stateOfMatter = Interfaces.Incompressible constrainedby Interfaces.StateOfMatter
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

     parameter Modelica.Units.SI.MolarFlowRate n_flow_assert(max=0) = -dropOfCommons.n_flow_reg "Assertion threshold for negative molar flows"
        annotation(Dialog(tab="Advanced"));

    Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

    Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

     Real r;

    protected
      outer Chemical.DropOfCommons dropOfCommons;

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

      Real r_in,n_flow_in,u_in,h_in;
      Real r_out,n_flow_out,u_out,h_out;

    equation
      assert(n_flow_in > n_flow_assert, "Negative massflow at Volume inlet", dropOfCommons.assertionLevel);
      assert(-n_flow_out > n_flow_assert, "Positive massflow at Volume outlet", dropOfCommons.assertionLevel);
      assert(x > 0, "Molar fraction must be positive");


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
     u_out = stateOfMatter.chemicalPotentialPure(
       substanceData,
       temperature,
       pressure,
       electricPotential,
       moleFractionBasedIonicStrength)
       + Modelica.Constants.R*temperature*log(a)
       + z*Modelica.Constants.F*electricPotential;

     h_out = molarEnthalpy;
     der(n_flow_out)*L = r_out;
     der(n_flow_in)*L = r_in - r;

     u_out = u_in + r;

     if not useInlet then
       n_flow_in = 0;
       h_in = 0;
       u_in = u_out;
     end if;
     if not useOutlet then
       n_flow_out = 0;
     end if;

     annotation (
       Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialSubstance;

    partial model PartialSubstanceInSolution "Substance properties for components, where the substance is connected with the solution"

      Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

      extends Boundaries.Internal.PartialSubstance;

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

    partial model PartialSubstanceInSolutionWithAdditionalPorts "Substance properties for components, where the substance is connected with the solution"

      extends Boundaries.Internal.PartialSubstanceInSolution;

    Modelica.Units.SI.MolarFlowRate q
      "Molar flow rate of the substance into the component";

      Interfaces.SubstanceMassPort_a port_m "Substance mass fraction port" annotation (Placement(transformation(extent={{92,-110},{112,-90}})));

      Interfaces.SubstanceMolarityPort_a port_c annotation (Placement(transformation(extent={{90,90},{110,110}})));

    equation
      //molar mass flow
      q=(n_flow_in + n_flow_out + port_c.q + port_m.m_flow/stateOfMatter.molarMassOfBaseMolecule(substanceData));

      //substance mass fraction
      port_m.x_mass = solution.mj/solution.m;
      port_c.c = solution.nj/solution.V;

    end PartialSubstanceInSolutionWithAdditionalPorts;
  annotation (Documentation(info="<html>
<p>This package contains all internal functions, partials, and other (e.g. experimental) models for the Boundaries package.</p>
</html>"));
  end Internal;
  annotation (Documentation(revisions="<html>
<p><img src=\"modelica:/ThermofluidStream/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>

</html>", info="<html>
<p>The boundaries are Sorces and Sinks, as well as Volumes, that are conceptually a source and a sink with extra equations and act as loop breakers in closes cycles, and therefore are also boundaries.</p>
</html>"));
end Boundaries;
