within Chemical;
package Boundaries "This package contains boundary models for the stream."
extends Modelica.Icons.SourcesPackage;

  model Substance "Substance in solution"
    extends Icons.Substance;
    extends Chemical.Boundaries.Internal.PartialSubstanceInSolution;

    Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
      "Molar concentration of particles";

     parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

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

    Real logn(stateSelect=StateSelect.prefer, start=log(m_start/molarMassOfBaseMolecule), min=0)
    "Natural logarithm of the amount of base molecules in solution";

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  initial equation

    amountOfBaseMolecules = m_start/molarMassOfBaseMolecule;

  equation
   substanceDataVar = substanceData;
   n_flow = n_flow_out;

    if stateOfMatter.selfClustering(substanceData) then

      //Liquid cluster theory - equilibrium:
      //substance.x[i] = substance.x*(K*substance.x)^i .. mole fraction of cluster composed with i base molecules
      //amountOfParticles/solution.n = substance.x/(1-K*substance.x);                //sum(substance.x[i])
      //amountOfBaseMolecules/solution.n = substance.x/((1-K*substance.x)^2);            //sum(i*substance.x[i])
      //amountOfHydrogenBonds/solution.n = substance.x*substance.x*K/((1-K*substance.x)^2);   //sum((i-1)*substance.x[i])

      amountOfParticles*(1 - SelfClustering_K*substance.x) = amountOfFreeMolecule;

      //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*substance.x)) = amountOfParticles":
      substance.x = ((2*SelfClustering_K+solution.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solution.n/amountOfBaseMolecules)+(solution.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

      amountOfBonds = amountOfBaseMolecules*substance.x*SelfClustering_K;

      //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
      //TODO: more precise calculation of other properties

     //der(enthalpy) = solution.dH + n_flow*actualStream(port_a.h_outflow);
     //enthalpy = substance.h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
      solution.dH =if (EnthalpyNotUsed) then 0 else der(substance.h)*
        amountOfBaseMolecules + n_flow*substance.h - n_flow_out*state_out.h  + (
        if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
        substanceData)*der(amountOfBonds) else 0)
                      "heat transfer from other substances in solution [J/s]";

      solution.Gj =amountOfBaseMolecules*state_out.u + amountOfBonds*SelfClustering_dG
                      "Gibbs energy of the substance";

    else

      amountOfParticles = amountOfFreeMolecule;
      amountOfBaseMolecules = amountOfFreeMolecule;
      amountOfBonds = 0;

      //der(enthalpy) = solution.dH + n_flow*actualStream(port_a.h_outflow);
      //enthalpy = substance.h*amountOfBaseMolecules;
      solution.dH =
        if (EnthalpyNotUsed) then  0
        else    der(substance.h)*amountOfBaseMolecules + n_flow*substance.h
                - n_flow_out*state_out.h
                "heat transfer from other substances in solution [J/s]";

      solution.Gj = amountOfBaseMolecules*state_out.u "Gibbs energy of the substance [J]";

    end if;

    //The main accumulation equation is "der(amountOfBaseMolecules)=n_flow"
    // However, the numerical solvers can handle it in form of log(n) much better. :-)
    der(logn) = (n_flow/amountOfBaseMolecules) "accumulation of amountOfBaseMolecules=exp(logn) [mol]";
    //der(amountOfBaseMolecules) = n_flow;
    amountOfBaseMolecules = exp(logn);

    substance.x = amountOfFreeMolecule/solution.n "mole fraction [mol/mol]";

    c = amountOfParticles/solution.V "concentration [mol/m3]";

    //solution flows
    solution.i = Modelica.Constants.F*substance.z*n_flow +
        Modelica.Constants.F*der(substance.z)*amountOfBaseMolecules "change of sunstance charge [A]";
    solution.dV = substance.Vm*n_flow + der(substance.Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

    //extensive properties
    solution.nj = amountOfParticles;
    solution.mj = amountOfBaseMolecules*molarMassOfBaseMolecule;
    solution.Vj = amountOfBaseMolecules*substance.Vm;
    solution.Qj = Modelica.Constants.F*amountOfBaseMolecules*substance.z;
    solution.Ij = (1/2)*(amountOfBaseMolecules*substance.z^2);

       annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={128,0,255},
            textString="%name")}), Documentation(revisions="<html>
<p>2009-2025 by Marek Matejak, Ph.D. </p>
</html>",   info="<html>
<h4>n = substance.x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and substance.x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
<p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
<p>substance.x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
<p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
<p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of substance.x will be the same as molality.</p>
<p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of substance.x will be the same as molarity.</p>
<p><br><br>Definition of electro-chemical potential:</p>
<h4>u = u&deg; + R*T*ln(gamma*substance.x) + substance.z*F*v</h4>
<h4>u&deg; = DfG = DfH - T * DfS</h4>
<p>where</p>
<p>substance.x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. relative eletric potential of the solution</p>
<p>substance.z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
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

  model Accumulation "Substance in solution with inlet"
    extends Icons.Substance;
    extends Chemical.Boundaries.Internal.PartialSubstanceInSolution;

    import Chemical.Utilities.Types.InitializationMethods;

    parameter InitializationMethods initAmount = Chemical.Utilities.Types.InitializationMethods.state "Initialization method for amount of substance";
     // annotation(Dialog(tab= "Initialization", group="Molar flow"));

    Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
      "Molar concentration of particles";

    Chemical.Interfaces.Inlet inlet(
      redeclare package stateOfMatter=stateOfMatter,
      r=r_in,
      n_flow=n_flow_in)
      "The substance entering"
      annotation (Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

    parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=1
    "Initial amount of substance base molecules"
      annotation ( Dialog(group="Initialization", enable=(initAmount == InitializationMethods.state)));

    parameter Modelica.Units.SI.MolarFlowRate change_start=1
    "Initial change of substance base molecules"
      annotation ( Dialog(group="Initialization", enable=(initAmount == InitializationMethods.derivative)));

    Modelica.Units.SI.Mass mass=amountOfBaseMolecules*
        molarMassOfBaseMolecule "Mass";

    parameter Boolean calculateClusteringHeat = true "Only for self clustering substances"
        annotation(Evaluate=true, choices(checkBox=true), Dialog(tab = "Clustering", enable = stateOfMatter.selfClustering(substanceDataVar)));

  protected
    Modelica.Units.SI.Mass m_start=amountOfSubstance_start*molarMassOfBaseMolecule;

    Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceDataVar);

    Modelica.Units.SI.AmountOfSubstance amountOfBaseMolecules(start=
         amountOfSubstance_start)
      "Amount of base molecules inside all clusters in compartment";

    Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule
      "Amount of free molecules not included inside any clusters in compartment";

    Modelica.Units.SI.AmountOfSubstance amountOfParticles
      "Amount of particles/clusters in compartment";

    Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
        Modelica.Constants.R*solution.T))
      "Dissociation constant of hydrogen bond between base molecules";

    Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
        stateOfMatter.selfClusteringBondEnthalpy(substanceDataVar)
      - solution.T * stateOfMatter.selfClusteringBondEntropy(substanceDataVar)
      "Gibbs energy of hydrogen bond between H2O molecules";

    Modelica.Units.SI.AmountOfSubstance amountOfBonds
      "Amount of hydrogen bonds between molecules in compartment";

    Real logn(stateSelect=StateSelect.prefer, start=log(amountOfSubstance_start), min=0)
    "Natural logarithm of the amount of base molecules in solution";

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

    Real r, r_in, n_flow_in;

    Chemical.Interfaces.SubstanceState state_in;


  initial equation
    if initAmount == InitializationMethods.steadyState then
      r=0;
    elseif initAmount == InitializationMethods.state then
      amountOfBaseMolecules = amountOfSubstance_start;
    elseif initAmount == InitializationMethods.derivative then
      n_flow = change_start;
    end if;

  equation
    substanceDataVar = inlet.definition;
    state_in = inlet.state;

    n_flow = n_flow_in + n_flow_out;

    der(n_flow_in)*L = r_in - r;
    state_out.u = state_in.u + r;

  if stateOfMatter.selfClustering(substanceDataVar) then

      //Liquid cluster theory - equilibrium:
      //substance.x[i] = substance.x*(K*substance.x)^i .. mole fraction of cluster composed with i base molecules
      //amountOfParticles/solution.n = substance.x/(1-K*substance.x);                //sum(substance.x[i])
      //amountOfBaseMolecules/solution.n = substance.x/((1-K*substance.x)^2);            //sum(i*substance.x[i])
      //amountOfHydrogenBonds/solution.n = substance.x*substance.x*K/((1-K*substance.x)^2);   //sum((i-1)*substance.x[i])

      amountOfParticles*(1 - SelfClustering_K*substance.x) = amountOfFreeMolecule;

      //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*substance.x)) = amountOfParticles":
      substance.x = ((2*SelfClustering_K+solution.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solution.n/amountOfBaseMolecules)+(solution.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

      amountOfBonds = amountOfBaseMolecules*substance.x*SelfClustering_K;

      //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
      //TODO: more precise calculation of other properties

     //der(enthalpy) = solution.dH + n_flow*actualStream(port_a.h_outflow);
     //enthalpy = substance.h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
      solution.dH =if (EnthalpyNotUsed) then 0 else der(substance.h)*
        amountOfBaseMolecules + n_flow*substance.h - n_flow_out*state_out.h - n_flow_in*state_in.h + (
        if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
        substanceDataVar)*der(amountOfBonds) else 0)
                      "heat transfer from other substances in solution [J/s]";

      solution.Gj =amountOfBaseMolecules*state_out.u + amountOfBonds*SelfClustering_dG
                      "Gibbs energy of the substance";

    else

      amountOfParticles = amountOfFreeMolecule;
      amountOfBaseMolecules = amountOfFreeMolecule;
      amountOfBonds = 0;

      //der(enthalpy) = solution.dH + n_flow*actualStream(port_a.h_outflow);
      //enthalpy = substance.h*amountOfBaseMolecules;
      solution.dH =
        if (EnthalpyNotUsed) then  0
        else    der(substance.h)*amountOfBaseMolecules + n_flow*substance.h
                - n_flow_out*state_out.h - n_flow_in*state_in.h
                "heat transfer from other substances in solution [J/s]";

      solution.Gj = amountOfBaseMolecules*state_out.u "Gibbs energy of the substance [J]";

    end if;

    //The main accumulation equation is "der(amountOfBaseMolecules)=n_flow"
    // However, the numerical solvers can handle it in form of log(n) much better. :-)
    der(logn) = (n_flow/amountOfBaseMolecules) "accumulation of amountOfBaseMolecules=exp(logn) [mol]";
    //der(amountOfBaseMolecules) = n_flow;
    amountOfBaseMolecules = exp(logn);

    substance.x = amountOfFreeMolecule/solution.n "mole fraction [mol/mol]";

    c = amountOfParticles/solution.V "concentration [mol/m3]";

    //solution flows
    solution.i = Modelica.Constants.F*substance.z*n_flow +
        Modelica.Constants.F*der(substance.z)*amountOfBaseMolecules "change of sunstance charge [A]";
    solution.dV = substance.Vm*n_flow + der(substance.Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

    //extensive properties
    solution.nj = amountOfParticles;
    solution.mj = amountOfBaseMolecules*molarMassOfBaseMolecule;
    solution.Vj = amountOfBaseMolecules*substance.Vm;
    solution.Qj = Modelica.Constants.F*amountOfBaseMolecules*substance.z;
    solution.Ij = (1/2)*(amountOfBaseMolecules*substance.z^2);

       annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={128,0,255},
            textString="%name")}), Documentation(revisions="<html>
<p>2009-2025 by Marek Matejak, Ph.D. </p>
</html>",   info="<html>
<h4>n = substance.x &middot; n(solution) = &int; MolarFlow</h4>
<p>where n is amount of the substance and substance.x is mole fraction.</p>
<p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
<p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
<p>substance.x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
<p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
<p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of substance.x will be the same as molality.</p>
<p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of substance.x will be the same as molarity.</p>
<p><br><br>Definition of electro-chemical potential:</p>
<h4>u = u&deg; + R*T*ln(gamma*substance.x) + substance.z*F*v</h4>
<h4>u&deg; = DfG = DfH - T * DfS</h4>
<p>where</p>
<p>substance.x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. relative eletric potential of the solution</p>
<p>substance.z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg; .. chemical potential of pure substance</p>
<p>DfG .. free Gibbs energy of formation of the substance</p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
  end Accumulation;

  model ElectronSource "Electron transfer from the solution to electric circuit"
    extends Icons.ElectronTransfer;

    Chemical.Interfaces.Outlet outlet(
      r=r_out,
      n_flow=n_flow,
      state(u=u, h=h),
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

    Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

    parameter Interfaces.Incompressible.SubstanceDataParameters substanceData = Chemical.Substances.Electrone_solid() "Definition of the substance";

    Real r_out, h;

    Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

    Modelica.Units.SI.ElectricPotential electricPotential "Electric potential of the solution";

    Modelica.Units.SI.Temperature temperature "Temperature of the solution";

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L;
    parameter Modelica.Units.SI.ChemicalPotential u_0=0 "Initial electro-chemical potential";

    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  initial equation
    u = u_0;
  equation

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
  end ElectronSource;

  model ElectronSink "Electron transfer to an electric circuit"
    extends Icons.ElectronTransfer;

    Chemical.Interfaces.Inlet inlet(
      n_flow=n_flow)
       "Chemical electron inlet" annotation (Placement(transformation(extent={{110,-10},{90,10}}), iconTransformation(extent={{110,-10},{90,10}})));

    Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
        Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
            extent={{-10,88},{10,108}})));

    Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

   parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L;

    Interfaces.Incompressible.SubstanceData substanceData;
    // == Chemical.Substances.Electrone_solid() "Definition of the substance";

    Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

    Modelica.Units.SI.ElectricPotential electricPotential "Electric potential of the solution";

    Modelica.Units.SI.Temperature temperature "Temperature of the solution";


    Modelica.Units.SI.ChemicalPotential u,r;

  protected
    outer Chemical.DropOfCommons dropOfCommons;

  initial equation
    r=0;
  equation

    substanceData=inlet.definition;


     der(inlet.n_flow)*L = inlet.r - r;
     u = r + inlet.state.u;


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

  model ExternalIdealGas "Ideal gas substance with defined partial pressure"
    extends Internal.PartialSubstanceInSolution(redeclare package stateOfMatter
        = Interfaces.IdealGas);
    extends Internal.PartialSolutionSensor;

     parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Boolean usePartialPressureInput = false
    "=true, if fixed partial pressure is from input instead of parameter"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Modelica.Units.SI.Pressure PartialPressure=0
      "Fixed partial pressure if usePartialPressureInput=false" annotation (
       HideResult=true, Dialog(enable=not usePartialPressureInput));

    Modelica.Blocks.Interfaces.RealInput partialPressure(start=
          PartialPressure, final unit="Pa")=p if usePartialPressureInput
    "Partial pressure of gas = total pressure * gas fraction"
      annotation (HideResult=true,Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.Pressure p "Current partial pressure";

  equation
   substanceDataVar = substanceData;
   n_flow = n_flow_out;

    if not usePartialPressureInput then
      p=PartialPressure;
    end if;

    //mole fraction
    substance.x = p / solution.p;

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
  end ExternalIdealGas;

  model ExternalPureSubstance "Constant source of pure substance"
    extends Internal.PartialSubstanceInSolution;
    extends Internal.PartialSolutionSensor;

    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

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
     substanceDataVar = substanceData;
     n_flow = n_flow_out;

     if stateOfMatter.selfClustering(substanceData) then

      //Liquid cluster theory - equilibrium:
      //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules

      //sum(x[i]) = x/(1-K*x) = amountOfParticles/amountOfParticles = 1;
      substance.x = 1/(1+SelfClustering_K) "mole fraction of free base molecule";
    else
      substance.x = 1 "pure substance is composed only with free base molecules";
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
  end ExternalPureSubstance;

   model ExternalMolality "Constant source of substance molality"
    extends Internal.PartialSubstanceInSolution;
    extends Internal.PartialSolutionSensor;

      parameter stateOfMatter.SubstanceDataParameters substanceData
    "Definition of the substance"
       annotation (choicesAllMatching = true);

     parameter Modelica.Units.SI.Molality Molality = 1e-8
    "Fixed molality of the substance if useMolalityInput=false"
      annotation (HideResult=true, Dialog(enable=not useMolalityInput));

      parameter Boolean useMolalityInput = false
    "Is amount of substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput molalityInput(start=Molality,final unit="mol/kg")=n/KG
      if useMolalityInput
      annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

   equation
     substanceDataVar = substanceData;
     n_flow = n_flow_out;


     if not useMolalityInput then
       n=Molality*solution.m;
     end if;

    substance.x = n/solution.n;

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
     extends Internal.PartialSubstanceInSolution;
     extends Internal.PartialSolutionSensor;

     parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);


     parameter Modelica.Units.SI.Concentration MolarConcentration = 1e-8
    "Fixed molarity of the substance if useMolarityInput=false"
      annotation (HideResult=true, Dialog(enable=not useMolarityInput));

    parameter Modelica.Units.SI.AmountOfSubstance AmountOfSolutionIn1L=55.508
      "Amount of all particles in the solution one liter of solvent";

      parameter Boolean useMolarityInput = false
    "Is amount of substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput molarConcentrationInput(start=MolarConcentration,final unit="mol/m3", displayUnit="mol/l")=n/L
      if useMolarityInput
      annotation (HideResult=true, Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Units.SI.AmountOfSubstance n "Current amount of the substance";

  equation
     substanceDataVar = substanceData;
     n_flow = n_flow_out;

     if not useMolarityInput then
       n=MolarConcentration*solution.V;
     end if;

    substance.x = n/solution.n;

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
       extends Internal.PartialSubstanceInSolution;
       extends Internal.PartialSolutionSensor;

    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Modelica.Units.SI.MoleFraction MoleFraction=1e-8
      "Fixed mole fraction of the substance if useMoleFractionInput=false"
      annotation (HideResult=true, Dialog(enable=not useMoleFractionInput));

      parameter Boolean useMoleFractionInput = false
    "Is mole fraction of the substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    Modelica.Blocks.Interfaces.RealInput moleFractionInput(
      final unit="mol/mol",
      start=MoleFraction)=x
      if useMoleFractionInput annotation (HideResult=true, Placement(transformation(
            extent={{-120,-20},{-80,20}})));

  equation
     substanceDataVar = substanceData;
     n_flow = n_flow_out;

     if not useMoleFractionInput then
       substance.x=MoleFraction;
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

  model ExternalChemicalPotentialSource "Constant source of electro-chemical potential"

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

    parameter Interfaces.SolutionStateParameters solutionState;
    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Boolean potentialFromInput = false "Use input connector for electro-chemical potential?";
    parameter Boolean enthalpyFromInput = false "Use input connector for molar enthalpy";
    parameter Boolean temperatureFromInput = false "Use input connector for temperature";

    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "Electro-chemical potential set value" annotation (Dialog(enable=not potentialFromInput));
    parameter Modelica.Units.SI.MolarEnthalpy h0_par=0 "Molar enthalpy set value"
      annotation (Dialog(enable=not enthalpyFromInput));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance" annotation (Dialog(tab="Advanced"));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J.mol-1") if potentialFromInput "Electro-chemical potential input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,40},{0,80}}), iconTransformation(extent={{-40,40},{0,80}})));
    Modelica.Blocks.Interfaces.RealInput h0_var(unit = "J/mol")
      if enthalpyFromInput "Enthalpy input connector [J/mol]"
      annotation (Placement(transformation(extent={{-40,-40},{0,0}}), iconTransformation(extent={{-42,-20},
              {-2,20}})));

    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(extent={{80,-20},{120,20}})));

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Blocks.Interfaces.RealInput u0(unit="J/mol") "Internal electro-chemical potential connector";
    Modelica.Blocks.Interfaces.RealInput h0(unit = "J/mol") "Internal enthalpy connector";

  equation
     outlet.definition = substanceData;
     outlet.solution = solutionState;

     connect(u0_var, u0);
     if not potentialFromInput then
       u0 = u0_par;
     end if;

     connect(h0_var, h0);
     if not enthalpyFromInput then
       h0 = h0_par;
     end if;


    L*der(outlet.n_flow) = outlet.r - 0;
    outlet.state.u = u0;
    outlet.state.h = h0;

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
  end ExternalChemicalPotentialSource;

  model ExternalChemicalPotentialSink "Boundary model of sink"

    parameter Boolean potentialFromInput = false "If true electro-chemical potential comes from real input";


    parameter Modelica.Units.SI.ChemicalPotential u0_par=0 "Electro-chemical potential setpoint of Sink" annotation (Dialog(enable=not potentialFromInput));
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of electro-chemical potential" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Blocks.Interfaces.RealInput u0_var(unit="J/mol")
      if potentialFromInput "Potential setpoint [J/mol]"
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

  //  Chemical.Utilities.Units.URT uRT=inlet.uRT;

  equation

    connect(u0_var, u0);
    if not potentialFromInput then
      u0 = u0_par;
    end if;


    der(inlet.n_flow)*L = inlet.r - r;
    r + inlet.state.u = u0;

  //  inlet.state.u0RT = uRT - log(x_par)/(Modelica.Constants.R*T0);

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
  end ExternalChemicalPotentialSink;

  model ExternalPartialPressureSink "Boundary model of sink"

    replaceable package stateOfMatterIn = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter "Substance model of inlet"
      annotation ( choices(
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    parameter Boolean partialPressureFromInput = false "If true partial pressure comes from real input";

    replaceable package stateOfMatter = Interfaces.IdealGas constrainedby
      Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
      annotation (choices(
       choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatter =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

     parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Modelica.Units.SI.Pressure p=1 "Partial pressure of substance - setpoint of Sink";
    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of electro-chemical potential" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter=stateOfMatterIn) annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));

    Modelica.Blocks.Interfaces.RealInput p_input(unit="Pa")=_p
      if partialPressureFromInput "Partial pressure of substance"
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

    Modelica.Units.SI.Pressure _p;

    Modelica.Units.SI.ChemicalPotential u,r,u_pure;

  equation

    if not partialPressureFromInput then
      _p = p;
    end if;

    der(inlet.n_flow)*L = inlet.r - r;
    r + inlet.state.u = u;

    u = stateOfMatter.chemicalPotentialPure(
     substanceData,
     inlet.solution.T,
     inlet.solution.p,
     inlet.solution.v,
     inlet.solution.I)
     + (Modelica.Constants.R*inlet.solution.T)*log(_p/inlet.solution.p)
     + inlet.definition.z*Modelica.Constants.F*inlet.solution.v;


    u_pure=stateOfMatter.chemicalPotentialPure(
     substanceData,
     inlet.solution.T,
     inlet.solution.p,
     inlet.solution.v,
     inlet.solution.I);

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
  end ExternalPartialPressureSink;

  model SubstanceInflow "Molar pump of substance to system"
    extends Interfaces.ConditionalSubstanceFlow;

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

    parameter Interfaces.SolutionStateParameters solutionState;
    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);


    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    Interfaces.Outlet outlet "Outflow" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_start;

  equation
    outlet.definition = substanceData;
    outlet.solution = solutionState;

    outlet.n_flow = q;

    TC * der(u) = outlet.r;
    outlet.state.u = u;
    outlet.state.h = stateOfMatter.molarEnthalpy(substanceData,solutionState.T);
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

    Interfaces.Inlet inlet "Inflow" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  equation
    inlet.n_flow = q;

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

  model TerminalSource "Source that imposes n_flow = 0"

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

    parameter Interfaces.SolutionStateParameters solutionState;
    parameter stateOfMatter.SubstanceDataParameters substanceData
   "Definition of the substance"
      annotation (choicesAllMatching = true);

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Modelica.Units.SI.MolarEnthalpy h=0 "Source enthalpy";
    parameter Modelica.Units.SI.ChemicalPotential u_start=0 "Initial electro-chemical potential";

    Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

  protected
    Modelica.Units.SI.ChemicalPotential u(stateSelect=StateSelect.prefer);

  initial equation
    u = u_start;

  equation
    outlet.definition = substanceData;
    outlet.solution = solutionState;

    outlet.n_flow = 0;

    TC * der(u) = outlet.r;
    outlet.state.u = u;
    outlet.state.h = h;

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

    Chemical.Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

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

  model OutStream "Outflow of whole solution"
    extends Boundaries.Internal.ConditionalSolutionFlow;
    extends Chemical.Interfaces.PartialSubstanceInlet;


  equation

    assert(volumeFlow>=-Modelica.Constants.eps, "Clearance can not be negative!");

    inlet.n_flow = c * volumeFlow;

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
  end OutStream;

  model Clearance "ClearanceFlow of whole solution"
    extends Chemical.Interfaces.PartialSubstanceInlet;

    parameter Modelica.Units.SI.VolumeFlowRate Clearance
    "Physiological clearance of the substance";

  equation

    assert(Clearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

    inlet.n_flow = c * Clearance;

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
    extends Chemical.Interfaces.PartialSubstanceInlet;

    parameter Modelica.Units.SI.Time HalfTime
    "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

    Modelica.Units.SI.AmountOfSubstance n;

  equation

    n = x*inlet.solution.n;
    inlet.n_flow = n * (Modelica.Math.log(2)/HalfTime);

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

  package Tests "Tests for the boundary package"
    extends Modelica.Icons.ExamplesPackage;

    model SourceSink "Test for source and sink model"
      extends Modelica.Icons.Example;

      inner Chemical.DropOfCommons dropOfCommons annotation (Placement(transformation(extent={{58,-74},{78,-54}})));
      ExternalChemicalPotentialSource source(
        u0_par=-200000, L=0,
        outlet(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-32,10},{-12,30}})));
      ExternalChemicalPotentialSink sink(u0_par=-100000) annotation (Placement(transformation(extent={{14,10},{34,30}})));
      ExternalChemicalPotentialSource source1(
        u0_par=-200000,
        outlet(n_flow(start=0, fixed=true)))
        annotation (Placement(transformation(extent={{-32,-30},{-12,-10}})));
      ExternalChemicalPotentialSink sink1(L=0, u0_par=-100000) annotation (Placement(transformation(extent={{14,-30},{34,-10}})));
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

     stateOfMatter.BaseProperties substance(substanceDataVar=substanceDataVar,T=temperature,p=pressure,v=electricPotential,I=moleFractionBasedIonicStrength);


     Modelica.Units.SI.Temperature temperature
      "Temperature of the solution";

    Modelica.Units.SI.Pressure pressure "Pressure of the solution";

    Modelica.Units.SI.ElectricPotential electricPotential
      "Electric potential of the solution";

    Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength
      "Ionic strength of the solution";



     outer Modelica.Fluid.System system "System wide properties";

     parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L;



      Chemical.Interfaces.Outlet outlet(
        redeclare package stateOfMatter = stateOfMatter,
        r=r_out,
        n_flow=n_flow_out,
        state=state_out,
        solution=solutionState,
        definition=substanceDataVar
        /*uRT=uRT_out,
    u0RT=u0RT,
    h=h_out*/)     "The substance exiting" annotation (Placement(transformation(extent={{90,-10},{110,10}})));




     stateOfMatter.SubstanceData substanceDataVar;

     parameter Modelica.Units.SI.MolarFlowRate n_flow_assert(max=0) = -dropOfCommons.n_flow_reg "Assertion threshold for negative molar flows"
        annotation(Dialog(tab="Advanced"));

     Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    /*
Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

Modelica.Units.SI.ActivityOfSolute a
  "Activity of the substance (mole-fraction based)";

 //Modelica.Units.SI.ChemicalPotential r "Inertial electro-chemical potential";

 
Modelica.Units.SI.ActivityCoefficient gamma
  "Activity coefficient of the substance";

Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";



//Modelica.Units.SI.MolarMass molarMass "Molar mass of the substance";

Modelica.Units.SI.MolarEnthalpy molarEnthalpy
  "Molar enthalpy of the substance";

Modelica.Units.SI.MolarEnthalpy molarEnthalpyPure
  "Molar enthalpy of the pure substance";

Modelica.Units.SI.MolarEntropy molarEntropyPure
  "Molar entropy of the pure substance";

Modelica.Units.SI.ChemicalPotential u0
  "Chemical potential of the pure substance";

Modelica.Units.SI.ChemicalPotential uPure
  "Electro-Chemical potential of the pure substance";

Modelica.Units.SI.MolarVolume molarVolume
  "Molar volume of the substance";

Modelica.Units.SI.MolarMass molarMass
  "Molar mass of the substance";

Modelica.Units.SI.MolarVolume molarVolumePure
  "Molar volume of the pure substance";

Modelica.Units.SI.MolarVolume molarVolumeExcess
  "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

*/
     Real r_out,n_flow_out;

     Chemical.Interfaces.SubstanceState state_out;

     Chemical.Interfaces.SolutionState solutionState;


    equation
      //assert(n_flow_in > n_flow_assert, "Negative massflow at Volume inlet", dropOfCommons.assertionLevel);
      //assert(-n_flow_out > n_flow_assert, "Positive massflow at Volume outlet", dropOfCommons.assertionLevel);
    /*  assert(x > 0, "Molar fraction must be positive");



 //aliases
 gamma = stateOfMatter.activityCoefficient(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 z = stateOfMatter.chargeNumberOfIon(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
// molarMass = stateOfMatter.molarMass(substanceDataVar);

 molarMass = 1/stateOfMatter.specificAmountOfParticles(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 molarEnthalpy = stateOfMatter.molarEnthalpy(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 molarEntropyPure = stateOfMatter.molarEntropyPure(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 molarEnthalpyPure = uPure + temperature*molarEntropyPure;
 u0 = stateOfMatter.chemicalPotentialPure(
   substanceDataVar,
   temperature,
   pressure,
   electricPotential,
   moleFractionBasedIonicStrength);
 uPure = stateOfMatter.electroChemicalPotentialPure(
   substanceDataVar,
   temperature,
   pressure,
   electricPotential,
   moleFractionBasedIonicStrength);
 molarVolume = stateOfMatter.molarVolume(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 molarVolumePure = stateOfMatter.molarVolumePure(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 molarVolumeExcess = stateOfMatter.molarVolumeExcess(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);
 //  molarHeatCapacityCp = stateOfMatter.molarHeatCapacityCp(substanceDataVar,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);

 //activity of the substance
 a = gamma*x;
*/
     //electro-chemical potential of the substance in the solution
     state_out.u = substance.u;
     /*stateOfMatter.chemicalPotentialPure(
   substanceDataVar,
   temperature,
   pressure,
   electricPotential,
   moleFractionBasedIonicStrength)
   + (Modelica.Constants.R*temperature)*log(a)
   + z*Modelica.Constants.F*electricPotential;
*/
     state_out.h = substance.h;

     der(n_flow_out)*L = r_out;
    /*
 if not useInlet then
   n_flow_in = 0;
   state_in = state_out;
   solutionState=solutionStateIn;
   substanceDataVarIn=substanceDataVar;
 end if;

*/
     annotation (
       Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Ph.D.</p>
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

      solutionState.T=temperature "Temperature of the solution";
      solutionState.p=pressure "Pressure of the solution";
      solutionState.v=electricPotential "Electric potential in the solution";
      solutionState.n=amountOfSolution "Amount of the solution";
      solutionState.m=solution.m "Mass of the solution";
      solutionState.V=solution.V "Volume of the solution";
      solutionState.G=solution.G "Free Gibbs energy of the solution";
      solutionState.Q=solution.Q "Electric charge of the solution";
      solutionState.I=solution.I "Mole fraction based ionic strength of the solution";

      //aliases
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;
      amountOfSolution = solution.n;
      moleFractionBasedIonicStrength = solution.I;

    end PartialSubstanceInSolution;

    partial model PartialSolutionSensor "Base class for sensor based on substance and solution properties"

      Interfaces.SolutionPort solution "To connect substance with solution, where is pressented"
        annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

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
  annotation (Documentation(info="<html>
<p>This package contains all internal functions, partials, and other (e.g. experimental) models for the Boundaries package.</p>
</html>"));
  end Internal;

  model Buffer
  "Source of substance bounded to constant amount of buffer to reach linear dependence between concentration and electrochemical potential"
  /*  extends Icons.Buffer;
  extends Internal.PartialSubstanceInSolution(
                 a(start = a_start));
  extends  Chemical.Obsolete.Interfaces.ConditionalKinetics(KC=1);


parameter Modelica.Units.SI.MoleFraction a_start=1e-7
  "Initial value of mole fraction of the buffered substance";

parameter Modelica.Units.SI.AmountOfSubstance BufferValue=0.001
  "Fixed buffer value (slope between amount of buffered substance and -log10(activity)) if useBufferValueInput=false"
  annotation (HideResult=true, Dialog(enable=not useBufferValueInput));

     parameter Boolean useBufferValueInput = false
  "Is buffer value of the substance an input?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));



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

    der(nFreeBuffer) = -n_flow;
    // <- This is mathematically the same as two following lines. However, the differential solvers can handle the log10n much better. :-)
    //der(log10nFreeBuffer)=(InvLog_10)*(port_a.q/nFreeBuffer);
    //nFreeBuffer = 10^log10nFreeBuffer;

    xFreeBuffer = nFreeBuffer/solution.n;
   // port_a.q = (solution.n*KC)*(xFreeBuffer - xref);
    n_flow = KC*(log(xFreeBuffer) - log(xref)); //alternative kinetics
    xref = -log10(a)*(bufferValue/solution.n);

  //der(n_flow)*L = r_in - r;
  //uRT = uRT_in + r;


  //solution flows
  streamEnthalpy = molarEnthalpy;

  solution.dH = state_out.h*n_flow_out - der(molarEnthalpy)*nFreeBuffer;
  solution.i = Modelica.Constants.F * z * n_flow - Modelica.Constants.F*der(z)*nFreeBuffer;
  solution.dV = molarVolume * n_flow - der(molarVolume)*nFreeBuffer;

  //extensive properties
  solution.nj=0;
  solution.mj=-nFreeBuffer*stateOfMatter.molarMassOfBaseMolecule(substanceData);
  solution.Vj=-nFreeBuffer*molarVolume;
  solution.Gj=-nFreeBuffer*state_out.u;
  solution.Qj=-Modelica.Constants.F*nFreeBuffer*z;
  solution.Ij=-(1/2) * ( nFreeBuffer * z^2);
*/
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
  annotation (Documentation(revisions="<html>
<p><img src=\"modelica:/ThermofluidStream/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>

</html>", info="<html>
<p>The boundaries are Sorces and Sinks, as well as Volumes, that are conceptually a source and a sink with extra equations and act as loop breakers in closes cycles, and therefore are also boundaries.</p>
</html>"));
end Boundaries;
