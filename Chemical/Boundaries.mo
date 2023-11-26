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

     parameter Real L=dropOfCommons.L;

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
