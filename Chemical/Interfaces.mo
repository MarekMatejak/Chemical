within Chemical;
package Interfaces "Chemical interfaces"
  import Chemical;

  connector Fore "Undirected connector outputting the forward state"

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.OutputSubstanceState state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.InputSubstanceState state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.OutputSolutionState solution "State of solution";
    Chemical.Interfaces.OutputDefinition definition "Definition of substance";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={194,138,221},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid)}),
                              Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Fore is a output. </u>
</html>"));

  end Fore;

  connector Rear "Undirected connector outputting the rearward state"

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.InputSubstanceState state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.OutputSubstanceState state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.InputSolutionState solution "State of solution";
    Chemical.Interfaces.InputDefinition definition "Definition of substance";

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={194,138,221},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
                              Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Rear is a input.</u>
</html>"));
  end Rear;

  extends Modelica.Icons.InterfacesPackage;

  operator record Definition "Definition of a chemical substance or a chemical process"

    Chemical.Interfaces.DataRecord data "Data record of the substance or process";



    encapsulated operator 'constructor'
      import Definition=Chemical.Interfaces.Definition;
      import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
      import DataRecord=Chemical.Interfaces.DataRecord;
      import PhaseType=Chemical.Interfaces.Phase;
      constant Real R=1.380649e-23*6.02214076e23;
      constant Real T0=298.15 "Base temperature";
      constant Real p0=100000 "Base pressure";

      function fromDataRecord
        input DataRecord data "Mass based data record";
        output Definition result(data=data) "Molar based data record";
      algorithm
        annotation (Inline=true);
      end fromDataRecord;

      function fromFormationEnergies
        input Real MM=1 "Molar mass of the substance";
        input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
        input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (298.15, 1bar)";
        input PhaseType phase=PhaseType.Incompressible "State of matter";
        input Real Vm=if (phase == PhaseType.Gas) then (R*T0)/p0 else 0.001*MM "Molar volume of the pure substance at SATP conditions (298.15, 1bar) (default 1L/mol)";
        input Real gamma=1 "Activity coefficient of the substance";
        output Definition result(
                  data=DataRecord(
                    MM=MM,
                    Hf=DfH,
                    H0=DfH - T0*Cp,
                    alow={0,0,Cp/R,0,0,0,0},
                    blow={(DfH-Cp*T0)/R,(((DfH-DfG)/T0)-Cp*log(T0))/R},
                    ahigh={0,0,Cp/R,0,0,0,0},
                    bhigh={(DfH-Cp*T0)/R,(((DfH-DfG)/T0)-Cp*log(T0))/R},
                    z=z,
                    phase=phase,
                    VmBase=Vm/(1+log(gamma)),
                    VmExcess=Vm*log(gamma)/(1+log(gamma))));
      algorithm
        annotation (Inline=true);
      end fromFormationEnergies;
    end 'constructor';

    encapsulated operator function '+'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;
      input Definition d1;
      input Definition d2;
      output Definition result " = d1 + d2";

    algorithm
      result :=Definition(
          data=DataRecord(
            MM=d1.data.MM + d2.data.MM,
            Hf=d1.data.Hf + d2.data.Hf,
            H0=d1.data.H0 + d2.data.H0,
            alow=d1.data.alow .+ d2.data.alow,
            blow=d1.data.blow .+ d2.data.blow,
            ahigh=d1.data.ahigh .+ d2.data.ahigh,
            bhigh=d1.data.bhigh .+ d2.data.bhigh,
            z=d1.data.z + d2.data.z,
            phase=d1.data.phase,
            VmBase=d1.data.VmBase + d2.data.VmBase,
            VmExcess=d1.data.VmExcess + d2.data.VmExcess));
          annotation (Inline=true);
    end '+';

    encapsulated operator '-'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;
     function negate
       input Definition d;
       output Definition result " = - d";
     algorithm
       result :=Definition(
          data=DataRecord(
            MM=-d.data.MM,
            Hf=-d.data.Hf,
            H0=-d.data.H0,
            alow=(-1) .* d.data.alow,
            blow=(-1) .* d.data.blow,
            ahigh=(-1) .* d.data.ahigh,
            bhigh=(-1) .* d.data.bhigh,
            z= -d.data.z,
            phase = d.data.phase,
            VmBase= -d.data.VmBase,
            VmExcess=d.data.VmExcess));

          annotation (Inline=true);
     end negate;

     function substract
      input Definition d1;
      input Definition d2;
      output Definition result " = d1 - d2";
     algorithm
      result :=Definition(
          data=DataRecord(
            MM=d1.data.MM - d2.data.MM,
            Hf=d1.data.Hf - d2.data.Hf,
            H0=d1.data.H0 - d2.data.H0,
            alow=d1.data.alow   .-  d2.data.alow,
            blow=d1.data.blow   .-  d2.data.blow,
            ahigh=d1.data.ahigh .-  d2.data.ahigh,
            bhigh=d1.data.bhigh .-  d2.data.bhigh,
            z=d1.data.z - d2.data.z,
            phase=d1.data.phase,
            VmBase=d1.data.VmBase - d2.data.VmBase,
            VmExcess=d1.data.VmExcess - d2.data.VmExcess));

          annotation (Inline=true);
     end substract;
    end '-';

    encapsulated operator '*'
      import Definition=Chemical.Interfaces.Definition;
      import DataRecord=Chemical.Interfaces.DataRecord;

    function scalar
      input Real n=1 "Stoichiometric coefficient";
      input Definition d;
      output Definition result " = n * d";
    algorithm
      result :=Definition(
          data=DataRecord(
            MM=n * d.data.MM,
            Hf=n * d.data.Hf,
            H0=n * d.data.H0,
            alow=n * d.data.alow,
            blow=n * d.data.blow,
            ahigh=n * d.data.ahigh,
            bhigh=n * d.data.bhigh,
            z=n * d.data.z,
            phase=d.data.phase,
            VmBase=n * d.data.VmBase,
            VmExcess=n * d.data.VmExcess));

          annotation (Inline=true);
    end scalar;

    function vector
      input Real[:] n "Stoichiometric coefficients";
      input Definition[:] d;
      output Definition result " = n * d";
    algorithm
      result :=Definition(
          data=DataRecord(
            MM=n * d.data.MM,
            Hf=n * d.data.Hf,
            H0=n * d.data.H0,
            alow= {sum({(n[i] * d[i].data.alow[j])  for i in 1:size(n,1)}) for j in 1:7},
            blow= {sum({(n[i] * d[i].data.blow[j])  for i in 1:size(n,1)}) for j in 1:2},
            ahigh={sum({(n[i] * d[i].data.ahigh[j]) for i in 1:size(n,1)}) for j in 1:7},
            bhigh={sum({(n[i] * d[i].data.bhigh[j]) for i in 1:size(n,1)}) for j in 1:2},
            z=n*d.data.z,
            phase=d[1].data.phase,
            VmBase=n*d.data.VmBase,
            VmExcess=n*d.data.VmExcess));
          annotation (Inline=true);
    end vector;
    end '*';

  end Definition;

  operator record SubstanceDefinition "Definition of a substance"
    //extends Chemical.Definition; //the inheritance of constructors is not supported in sufficient form in Modelica 3.6
    Chemical.Interfaces.DataRecord data "Data record of the base substance";

    Boolean SelfClustering "default=false If true then the base molecules are binding together into clusters";
    Modelica.Units.SI.MolarEnthalpy SelfClustering_dH "Enthalpy of bond between two base molecules of substance at 25degC, 1 bar";
    Modelica.Units.SI.MolarEntropy SelfClustering_dS "Entropy of bond between two base molecules of substance at 25degC, 1 bar";

    encapsulated operator 'constructor'
      import Definition=Chemical.Interfaces.SubstanceDefinition;
      import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
      import DataRecord=Chemical.Interfaces.DataRecord;
      import PhaseType=Chemical.Interfaces.Phase;
      constant Real R=1.380649e-23*6.02214076e23;
      constant Real T0=298.15 "Base temperature";
      constant Real p0=100000 "Base pressure";

      function fromDataRecord
        input DataRecord data "Mass based data record";
        input Boolean SelfClustering = false;
        input Real SelfClustering_dH = 0;
        input Real SelfClustering_dS = 0;
        output Definition result(
                  data=data,
                  SelfClustering=SelfClustering,
                  SelfClustering_dH=SelfClustering_dH,
                  SelfClustering_dS=SelfClustering_dS) "Molar based data record";
      algorithm
        annotation (Inline=true);
      end fromDataRecord;

      function fromFormationEnergies
        input Real MM=1 "Molar mass of the substance";
        input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
        input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (298.15, 1bar)";
        input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (298.15, 1bar)";
        input PhaseType phase=PhaseType.Incompressible "State of matter";
        input Real Vm=if (phase == PhaseType.Gas) then (R*T0)/p0 else 0.001*MM "Molar volume of the pure substance at SATP conditions (298.15, 1bar) (default fron non-gaseous is to reach density 1kg/L)";
        input Real gamma=1 "Activity coefficient of the substance";
        input Boolean SelfClustering = false;
        input Real SelfClustering_dH = 0;
        input Real SelfClustering_dS = 0;
        output Definition result(
                  data=DataRecord(
                    MM=MM,
                    Hf=DfH,
                    H0=DfH - T0*Cp,
                    alow={0,0,Cp/R,0,0,0,0},
                    blow={(DfH-Cp*T0)/R,(((DfH-DfG)/T0)-Cp*log(T0))/R},
                    ahigh={0,0,Cp/R,0,0,0,0},
                    bhigh={(DfH-Cp*T0)/R,(((DfH-DfG)/T0)-Cp*log(T0))/R},
                    z=z,
                    phase=phase,
                    VmBase=Vm/(1+log(gamma)),
                    VmExcess=Vm*log(gamma)/(1+log(gamma))),
                  SelfClustering=SelfClustering,
                  SelfClustering_dH=SelfClustering_dH,
                  SelfClustering_dS=SelfClustering_dS);
      algorithm
        annotation (Inline=true);
      end fromFormationEnergies;
    end 'constructor';
  end SubstanceDefinition;

  type Phase                 = enumeration(
    Gas
      "Gaseous phase",
    Liquid
      "Liquid phase",
    Solid
      "Liquid phase",
    Aqueous
      "Dissolved in water",
    Incompressible
      "Incompressible liquid or solid phase")
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>
Phase of substance or solution. It is possible to include any phase - as default it behave as incompressible. 
To change its behavior it is necessary to modify Property functions.
</p>
</html>"));
  package Properties "Calculation of chemical properties from definition and solution state"

    constant Real R=1.380649e-23*6.02214076e23;
    constant Real T0=298.15 "Base temperature";
    constant Real p0=100000 "Base pressure";

    model BaseSubstanceProperties "Base properties of the substance"
      import Chemical;

      parameter Boolean FixedDefinition "definition==definitionParam";

      parameter Chemical.Interfaces.SubstanceDefinition definitionParam "used only if FixedDefinition to help initialization";

      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

      parameter Boolean SolutionObserverOnly = false "True if substance does not affect the solution";

      Interfaces.InputSubstanceDefinition definition "Definition of the substance";

      Interfaces.InputSolutionState solutionState "State of the solution";

      InputAmountOfSubstance amountOfBaseMolecules
        "Amount of base molecules inside all clusters in compartment";

      InputMolarFlowRate n_flow "Molar change of the substance";

      InputHeatFlowRate h_flow "Substance enthaply change";

      Modelica.Units.SI.MoleFraction x "Mole fraction of the base molecule of the substance";

      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
            "Molar concentration of particles";

      Modelica.Units.SI.MassConcentration M(displayUnit="mg/l")
            "Mass concentration";

      Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
            "Molality";

      Modelica.Units.SI.MassFraction X
            "Mass fraction";

      Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";

      Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";

      Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of the substance";

      Modelica.Units.SI.MolarEnthalpy hPure
      "Molar enthalpy of the pure substance";

      Modelica.Units.SI.MolarEntropy sPure
      "Molar entropy of the pure substance";

      Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

      Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

      Modelica.Units.SI.MolarVolume Vm
      "Molar volume of the substance";

      Modelica.Units.SI.MolarMass MM
      "Molar mass of the substance";

      Modelica.Units.SI.MolarVolume VmPure
      "Molar volume of the pure substance";

      Modelica.Units.SI.MolarVolume VmExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

      Modelica.Units.SI.SpecificVolume Vs "Specific volume";

      Modelica.Units.SI.Density rho "Density";

     // Local connector definition, used for equation balancing check

      connector InputMolarFlowRate = input Modelica.Units.SI.MolarFlowRate
        "Molar flow rate as input signal connector";
      connector InputHeatFlowRate = input Modelica.Units.SI.HeatFlowRate
        "Heat flow rate as input signal connector";

      connector InputSubstanceData = input Chemical.Interfaces.Definition
        "Substance definition as input signal connector";
      connector InputTemperature = input Modelica.Units.SI.Temperature
        "Temperature as input signal connector";
      connector InputPressure = input Modelica.Units.SI.Pressure
        "Pressure as input signal connector";
      connector InputElectricPotential = input
          Modelica.Units.SI.ElectricPotential
        "Electric potential as input signal connector";
      connector InputMoleFraction = input Modelica.Units.SI.MoleFraction
        "Mole fraction as input signal connector";
      connector InputAmountOfSubstance = input
          Modelica.Units.SI.AmountOfSubstance
        "Amount of substance as input signal connector";

      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*specificAmountOfFreeBaseMolecule(
                                       definitionParam,
                                       Chemical.Interfaces.SolutionStateParameters(
                                       T=system.T_ambient,
                                       p=system.p_ambient)))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*specificAmountOfParticles(
                                       definitionParam,
                                       Chemical.Interfaces.SolutionStateParameters(
                                       T=system.T_ambient,
                                       p=system.p_ambient)))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solutionState.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          selfClusteringBondEnthalpy(definition)
        - solutionState.T * selfClusteringBondEntropy(definition)
        "Gibbs energy of hydrogen bond between H2O molecules";

      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";

      outer Modelica.Fluid.System system "System wide properties";

      Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;

    equation
     // assert(x > 0, "Molar fraction must be positive");

     gamma = activityCoefficient(definition,solutionState);
     z = chargeNumberOfIon(definition,solutionState);

     MM = 1/specificAmountOfParticles(definition,solutionState);
     h = molarEnthalpy(definition,solutionState);
     sPure = molarEntropyPure(definition,solutionState);
     hPure = uPure + solutionState.T*sPure;
     u0 = chemicalPotentialPure(definition,solutionState);
     uPure = electroChemicalPotentialPure(definition,solutionState);
     Vm = molarVolume(definition,solutionState);
     VmPure =Chemical.Interfaces.Properties.molarVolumeBase(definition, solutionState);
     VmExcess = molarVolumeExcess(definition,solutionState);
     //  molarHeatCapacityCp = smolarHeatCapacityCp(definition,solutionState);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     u = chemicalPotentialPure(definition,solutionState)
       + (Modelica.Constants.R*solutionState.T)*log(a)
       + z*Modelica.Constants.F*solutionState.v;

      // during initialization it does not take value from parameter definition if not useInlet,
      // so instead of just "selfClustering(definition)" it must be written
      if (FixedDefinition and selfClustering(definitionParam)) or selfClustering(definition) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solutionState.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solutionState.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solutionState.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        x = ((2*SelfClustering_K+solutionState.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solutionState.n/amountOfBaseMolecules)+(solutionState.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties

       //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
       //enthalpy = h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        dH = der(h)*
          amountOfBaseMolecules + n_flow*h - h_flow +
          selfClusteringBondEnthalpy(definition)*der(amountOfBonds)
                        "heat transfer from other substances in solution [J/s]";

        Gj =amountOfBaseMolecules*u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
        //enthalpy = h*amountOfBaseMolecules;
        dH =  der(h)*amountOfBaseMolecules + n_flow*h - h_flow
                  "heat transfer from other substances in solution [J/s]";

        Gj = amountOfBaseMolecules*u "Gibbs energy of the substance [J]";

      end if;

      x = amountOfFreeMolecule/solutionState.n "mole fraction of base molecule [mol/mol]";

      c = amountOfParticles/solutionState.V "concentration [mol/m3]";

      M = (amountOfParticles/specificAmountOfParticles(definition,solutionState))/solutionState.V "mass concentration [kg/m3]";

      b = amountOfParticles/solutionState.m "molality [mol/kg]";

      X  = (amountOfParticles/specificAmountOfParticles(definition,solutionState))/solutionState.m "mass fraction [kg/kg]";

      Vs = specificVolume(definition,solutionState);

      rho = density(definition,solutionState);

      if SolutionObserverOnly then
        i = 0;
        dV = 0;
        nj = 0;
        mj = 0;
        Vj = 0;
        Qj = 0;
        Ij = 0;
      else
      //solution flows
        i = Modelica.Constants.F*z*n_flow +
          Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
        dV = Vm*n_flow + der(Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
        nj = amountOfParticles;
        mj = amountOfBaseMolecules*molarMassOfBaseMolecule(definition);
        Vj = amountOfBaseMolecules*Vm;
        Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
        Ij = (1/2)*(amountOfBaseMolecules*z^2);
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end BaseSubstanceProperties;

    model BaseProcessProperties "Base properties of the chemical process"
      import Chemical;

      Interfaces.InputDefinition definition "Definition of the process";

      Interfaces.InputSolutionState solutionState "State of the solution";

      Modelica.Units.SI.ChemicalPotential dG "Gibbs energy change during the process";

      Modelica.Units.SI.MolarEnthalpy dH "Molar enthalpy change during the process";

      Modelica.Units.SI.MolarEntropy dS "Molar entropy change during the process";

      Modelica.Units.SI.MolarHeatCapacity dCp "Molar heat capacity change during the process";

      Modelica.Units.SI.MolarVolume dVm "Molar volume change during the process";

      Modelica.Units.SI.MolarVolume dVmExcess "Molar volume excess change during the process";


      outer Modelica.Fluid.System system "System wide properties";


    equation
     assert(abs(definition.data.MM) < Modelica.Constants.eps, "Process should not change the mass");
     assert(abs(chargeNumberOfIon(definition,solutionState)) < Modelica.Constants.eps, "Process should not change the charge");

     dH = molarEnthalpy(definition,solutionState);
     dG = chemicalPotentialPure(definition,solutionState);
     dG = dH - solutionState.T*dS;//  dS = molarEntropy(definition,dG,solutionState);

     dCp = molarHeatCapacityCp(definition,solutionState);
     dVm = molarVolume(definition,solutionState);
     dVmExcess = molarVolumeExcess(definition,solutionState);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end BaseProcessProperties;

   function processData
     "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
     import Chemical;

        extends Modelica.Icons.Function;
    input Real K "Process dissociation constant (mole-fraction based) at 25°C,1bar";
    input Modelica.Units.SI.MolarEnergy dH=0 "Process molar enthalpy change at 25°C,1bar";
    input Modelica.Units.SI.MolarHeatCapacity dCp=0 "Process molar heat capacity change at 25°C,1bar";
    input Modelica.Units.SI.MolarVolume dVm=0 "Process molar volume change at 25°C,1bar";
     input Chemical.Interfaces.Phase phase=Chemical.Interfaces.Phase.Incompressible "State of matter";
     output Chemical.Interfaces.Definition processData "Data record of process changes";
   algorithm
       processData :=Chemical.Interfaces.Definition(
       MM=0,
       z=0,
       DfG=-Modelica.Constants.R*T0*log(K),
       DfH=dH,
       Cp=dCp,
       Vm=dVm,
       phase=phase,
       gamma=1);


   end processData;

   function firstProductDefinition
       "Return formation definition of the first product of chemical process"
     import Chemical;

        extends Modelica.Icons.Function;
    //input
    input Modelica.Units.SI.StoichiometricNumber s[:] "Stoichiometric reaction coefficient for substrates [nS]";
    input Modelica.Units.SI.StoichiometricNumber p[:] "Stoichiometric reaction coefficient for products [nP]";
     input Chemical.Interfaces.Definition process "Data record of process changes";
     input Chemical.Interfaces.Definition substrates[:] "Substrates definitions [nS]";
     input Chemical.Interfaces.SubstanceDefinition products[:] "Products definitions [nP]";
     output Chemical.Interfaces.Definition firstProductDefinition "Definition of the first product in process";
   protected
     Chemical.Interfaces.Definition pd[size(products,1)-1];
   algorithm
        for i in 1:size(products,1)-1 loop
          pd[i].data := products[i+1].data;
        end for;
        firstProductDefinition := (1/p[1]) * (s*substrates + process - p[2:end]*pd);
         /*Definition(
       MM = ((p[2:end]*products.data.MM) - (s*substrates.MolarWeight))/p[1],
       z = ((p[2:end]*products.z) - (s*substrates.z))/p[1],
       DfG = ((p[2:end]*products.DfG) - (s*substrates.DfG) - process.DfG)/p[1],
       DfH = ((p[2:end]*products.DfH) - (s*substrates.DfH) - process.DfH)/p[1],
       gamma = 1,
       Cp = ((p[2:end]*products.Cp) - (s*substrates.Cp) - process.Cp)/p[1],
       SelfClustering = 0,
       SelfClustering_dH = 0,
       SelfClustering_dS = 0,
       Vs = ((p[2:end]*products.Vs) - (s*substrates.Vs) - process.Vs)/p[1]);*/
   end firstProductDefinition;

   function activityCoefficient
    "Return activity coefficient of the substance in the solution"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Real activityCoefficient "Activity Coefficient";
   algorithm
       activityCoefficient := if definition.data.phase==Chemical.Interfaces.Phase.Gas then 1 else exp(definition.data.VmExcess/definition.data.VmBase);
   end activityCoefficient;

   function chargeNumberOfIon
    "Return charge number of the substance in the solution"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
      "Charge number of ion";
   algorithm
       chargeNumberOfIon := definition.data.z;
   end chargeNumberOfIon;

   function molarEnthalpyElectroneutral
    "Molar enthalpy of the substance in electroneutral solution"
     import Chemical;
      extends Modelica.Icons.Function;
      import Modelica.Math;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

      output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
      "Molar enthalpy";
   algorithm
       //Molar enthalpy:
       // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

       molarEnthalpyElectroneutral := smooth(0,(if solution.T < Tlimit then R*((-definition.data.alow[1] + solution.T*(definition.data.blow[1] +
       definition.data.alow[2]*Math.log(solution.T) + solution.T*(1.*definition.data.alow[3] + solution.T*(0.5*definition.data.alow[4] +
       solution.T*(1/3*definition.data.alow[5] + solution.T*(0.25*definition.data.alow[6] + 0.2*definition.data.alow[7]*solution.T))))))
       /solution.T) else R*((-definition.data.ahigh[1] + solution.T*(definition.data.bhigh[1] + definition.data.ahigh[2]*
       Math.log(solution.T) + solution.T*(1.*definition.data.ahigh[3] + solution.T*(0.5*definition.data.ahigh[4] + solution.T*(1/3*definition.data.ahigh[5] +
       solution.T*(0.25*definition.data.ahigh[6] + 0.2*definition.data.ahigh[7]*solution.T))))))/solution.T)));


   end molarEnthalpyElectroneutral;

   function molarEnthalpy
    "Molar enthalpy of the substance with electric potential dependence"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      output Modelica.Units.SI.MolarEnthalpy molarEnthalpy "Molar enthalpy";
   algorithm
      molarEnthalpy := molarEnthalpyElectroneutral(definition,solution) +
           Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
      annotation (Inline=true, smoothOrder=2);
   end molarEnthalpy;

   function molarEntropyPure
    "Molar entropy of the pure substance"
     import Chemical;
      extends Modelica.Icons.Function;
      import Modelica.Math;
     input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";

   algorithm
       //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
       //+ substanceData.Cv*log(solution.T/298.15);

       //Molar entropy shift:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = solution.T*dS (small amount of added heat energy)
       // - pressure shift: with constant molar volume at constant temperature Vm*dP = -solution.T*dS (small amount of work)
       molarEntropyPure :=(if solution.T < Tlimit then R*(definition.data.blow[2] - 0.5*definition.data.alow[1]/(solution.T*solution.T) - definition.data.alow[
        2]/solution.T + definition.data.alow[3]*Math.log(solution.T) + solution.T*(definition.data.alow[4] + solution.T*(0.5*definition.data.alow[5] + solution.T
        *(1/3*definition.data.alow[6] + 0.25*definition.data.alow[7]*solution.T)))) else R*(definition.data.bhigh[2] - 0.5*definition.data.ahigh[1]/(solution.T
        *solution.T) - definition.data.ahigh[2]/solution.T + definition.data.ahigh[3]*Math.log(solution.T) + solution.T*(definition.data.ahigh[4] + solution.T*
        (0.5*definition.data.ahigh[5] + solution.T*(1/3*definition.data.ahigh[6] + 0.25*definition.data.ahigh[7]*solution.T))))) + (if definition.data.phase
         == Chemical.Interfaces.Phase.Gas then -R*log(solution.p/100000) else -(Chemical.Interfaces.Properties.molarVolumeBase(definition, solution)/solution.T)
        *(solution.p - 100000));


       //For example at triple point of water should be solution.T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //As definition.data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
       //At solution.T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K

   end molarEntropyPure;

    function molarEntropy "Molar entropy of the substance in the solution"
      import Chemical;
          extends Modelica.Icons.Function;
      input Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance";
      input Chemical.Interfaces.Definition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";

      output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
    algorithm
        molarEntropy :=  (u - molarEnthalpy(definition,solution))/solution.T;
    end molarEntropy;

   function chemicalPotentialPure "Chemical potential of the pure substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
      "Base chemical potential";
   algorithm
       chemicalPotentialPure :=  molarEnthalpyElectroneutral(definition,solution) - solution.T*molarEntropyPure(definition,solution);
   end chemicalPotentialPure;

   function electroChemicalPotentialPure
    "Electro-chemical potential of the pure substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.ChemicalPotential
      electroChemicalPotentialPure "Base electro-chemical potential";
   algorithm
    electroChemicalPotentialPure := chemicalPotentialPure(
         definition,
         solution) + Modelica.Constants.F*chargeNumberOfIon(definition,solution)*solution.v;
   end electroChemicalPotentialPure;

   function molarVolumeBase "Base molar volume of the substance"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
   algorithm
     molarVolumePure :=(if definition.data.phase == Chemical.Interfaces.Phase.Gas then R*solution.T/solution.p else definition.data.VmBase);

   end molarVolumeBase;

   function molarVolumeExcess
    "Excess molar volume of the substance in the solution"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Excess molar volume of the substance in the solution";
   algorithm
     //zero if activityCoefficient==1
      molarVolumeExcess :=(if definition.data.phase == Chemical.Interfaces.Phase.Gas then 0 else definition.data.VmExcess);

      annotation (Inline=true, smoothOrder=2);
   end molarVolumeExcess;

   function molarVolume "Molar volume of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
   algorithm
    molarVolume :=Chemical.Interfaces.Properties.molarVolumeBase(definition, solution) + molarVolumeExcess(definition, solution);
      annotation (Inline=true, smoothOrder=2);
   end molarVolume;

   function molarHeatCapacityCp
    "Molar heat capacity at constant pressure"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
      "Molar heat capacity at constant pressure";
   algorithm
     molarHeatCapacityCp := smooth(0,if solution.T < Tlimit then R*(1/(solution.T*solution.T)*(definition.data.alow[1] + solution.T*(
       definition.data.alow[2] + solution.T*(1.*definition.data.alow[3] + solution.T*(definition.data.alow[4] + solution.T*(definition.data.alow[5] + solution.T
       *(definition.data.alow[6] + definition.data.alow[7]*solution.T))))))) else R*(1/(solution.T*solution.T)*(definition.data.ahigh[1]
        + solution.T*(definition.data.ahigh[2] + solution.T*(1.*definition.data.ahigh[3] + solution.T*(definition.data.ahigh[4] + solution.T*(definition.data.
       ahigh[5] + solution.T*(definition.data.ahigh[6] + definition.data.ahigh[7]*solution.T))))))));
   end molarHeatCapacityCp;

   function molarMassOfBaseMolecule
      "Molar mass of base molecule of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarMass molarMass "Molar mass";
   algorithm
     molarMass := definition.data.MM;
   end molarMassOfBaseMolecule;

   function selfClustering "returns true if substance molecules are joining together to clusters"
     import Chemical;
       extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Data record of substance";
          output Boolean selfClustering;
   algorithm
     selfClustering:=definition.SelfClustering;
   end selfClustering;

   function selfClusteringBondEnthalpy
    "Enthalpy of joining two base molecules of the substance together to cluster"
     import Chemical;
       extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Data record of substance";
    output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
   algorithm
     selfClusteringEnthalpy:=definition.SelfClustering_dH;
   end selfClusteringBondEnthalpy;

   function selfClusteringBondEntropy
    "Entropy of joining two base molecules of the substance together to cluster"
     import Chemical;
       extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Data record of substance";
    output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
   algorithm
     selfClusteringEntropy:=definition.SelfClustering_dS;
   end selfClusteringBondEntropy;

   function selfClusteringBondVolume
     import Chemical;
       extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Data record of substance";
    output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
   algorithm
     selfClusteringBondVolume:=0;
   end selfClusteringBondVolume;

   function selfClusteringBondHeatCapacityCp
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.Definition definition "Data record of substance";
    output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
   algorithm
     selfClusteringBondHeatCapacityCp:=0;
   end selfClusteringBondHeatCapacityCp;

    function specificAmountOfParticles
      "Amount of particles per mass of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.SubstanceDefinition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg")
        "Amount of substance particles per its mass";

    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
      Real amountOfBaseMolecules,amountOfFreeMolecule,amountOfParticles;
      Real x;
    algorithm
      if not selfClustering(definition) then
        specificAmountOfSubstance := 1/molarMassOfBaseMolecule(definition);
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));

        amountOfBaseMolecules:=mass/molarMassOfBaseMolecule(definition);
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfFreeMolecule := x*nSolution;

        amountOfParticles := amountOfFreeMolecule/(1 - SelfClustering_K*x);

        specificAmountOfSubstance := amountOfParticles/mass;

        //specificAmountOfSubstance := 1/((SelfClustering_K + 1)*definition.MolarWeight);
      end if;

      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfParticles;

    function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.SubstanceDefinition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";

    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K,amountOfBaseMolecules,x;
    algorithm
      if not selfClustering(definition) then
        specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(definition);
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));

        amountOfBaseMolecules:=mass/molarMassOfBaseMolecule(definition);
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        specificAmountOfFreeBaseMolecule := (x*nSolution)/mass;

      end if;
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

   function specificEnthalpy
     "Specific molar enthalpy of the substance with electric potential dependence"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
      "Specific enthalpy";

   protected
     Modelica.Units.SI.MolarEnergy SelfClustering_dG;
     Real SelfClustering_K;
   algorithm
     if selfClustering(definition) then
       SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
         selfClusteringBondEntropy(definition);
       SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
     end if;

     specificEnthalpy := (molarEnthalpy(definition,solution) +
       (if selfClustering(definition) then
       selfClusteringBondEnthalpy(definition)*SelfClustering_K/(
       SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

   /*algorithm 

  specificEnthalpy := molarEnthalpy(
    definition,
    solution)/
    molarMassOfBaseMolecule(definition);*/
   end specificEnthalpy;

   function specificVolume "Specific volume of the substance"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.SubstanceDefinition definition "Definition of substance";
     input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

   protected
     Modelica.Units.SI.MolarEnergy SelfClustering_dG;
     Real SelfClustering_K;
   algorithm
     if selfClustering(definition) then
       SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
         selfClusteringBondEntropy(definition);
       SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
     end if;

     specificVolume := (molarVolume(definition,solution) +
     (if selfClustering(definition) then
       selfClusteringBondVolume(definition)*SelfClustering_K/(
       SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

   end specificVolume;

    function specificHeatCapacityCp
    "Specific heat capacity at constant pressure"
      import Chemical;
      extends Modelica.Icons.Function;
      input Chemical.Interfaces.SubstanceDefinition definition "Definition of substance";
      input Interfaces.SolutionState solution "Chemical solution state";
    output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
      "Specific heat capacity at constant pressure";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(definition) then
        SelfClustering_dG := selfClusteringBondEnthalpy(definition) - solution.T*
          selfClusteringBondEntropy(definition);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*solution.T));
      end if;

      specificHeatCapacityCp := (molarHeatCapacityCp(
          definition,solution) + (if selfClustering(definition) then
        selfClusteringBondHeatCapacityCp(definition)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(definition);

        //TODO: + selfClusteringBondEnthalpy * der(K/(K + 1))/der(solution.T) .. if (selfClusteringBondHeatCapacityCp!=0)

    end specificHeatCapacityCp;

   function temperature
    "Temperature of the substance from its enthalpy"
     import Chemical;
      extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
     input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";

     input Modelica.Units.SI.Pressure p=100000 "Pressure";
     input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
     input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";

   protected
      function f_nonlinear "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
        input Interfaces.DataRecord data "Ideal gas data";
        input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";
      algorithm
        y := H_T(data=data, T=u) - h;
      end f_nonlinear;
   algorithm
      T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function f_nonlinear(data=definition.data, h=h), 200, 6000);

   /*protected 
  Modelica.Units.SI.SpecificEnthalpy baseSpecificEnthalpy;
algorithm 

  baseSpecificEnthalpy := specificEnthalpy(
      definition,
      298.15,
      p,
      v,
      I);

  T := 298.15 + (h - baseSpecificEnthalpy)/specificHeatCapacityCp(
  substanceData); */
   end temperature;

   function solution_temperature
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
     import Chemical;

       extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition[:] "Definition of substances";
    input Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of solution";
    input Modelica.Units.SI.MoleFraction x[:]
      "Mole fractions of substances";

    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";

   protected
     Chemical.Interfaces.Definition solutionDefinition=x*definition;

   algorithm
       T :=temperature(
           solutionDefinition,
           h,
           p,
           v,
           I);
   end solution_temperature;

   function density
        "Return density of the substance in the solution"
     import Chemical;
          extends Modelica.Icons.Function;
     input Chemical.Interfaces.Definition definition "Definition of substance";
    input Interfaces.SolutionState solution "Chemical solution state";

    output Modelica.Units.SI.Density density "Density";

   algorithm
     density :=if definition.data.phase == Chemical.Interfaces.Phase.Gas then (definition.data.MM*solution.p)/(R*solution.T) else definition.data.MM/(definition.data.VmBase+definition.data.VmExcess);

   end density;

   function H_T "Molar enthalpy of the substance in electroneutral solution"
      extends Modelica.Icons.Function;
      import Modelica.Math;
      input Interfaces.DataRecord data "Data record of the substance";
      input Modelica.Units.SI.Temperature T "Temperature";

      output Modelica.Units.SI.MolarEnthalpy H "Molar enthalpy";
   algorithm
       H :=
       smooth(0,(if T < Tlimit then R*((-data.alow[1] + T*(data.blow[1] +
       data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.alow[4] +
       T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
       /T) else R*((-data.ahigh[1] + T*(data.bhigh[1] + data.ahigh[2]*
       Math.log(T) + T*(1.*data.ahigh[3] + T*(0.5*data.ahigh[4] + T*(1/3*data.ahigh[5] +
       T*(0.25*data.ahigh[6] + 0.2*data.ahigh[7]*T))))))/T)) +
       data.H0);
       //(if   exclEnthForm then -data.Hf else 0.0) +
       //(if (refChoice == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) +
       //(if refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else 0.0))

   end H_T;
    annotation (Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Properties;

  constant Modelica.Units.SI.Temperature Tlimit=1000 "Temperature limit between low and high data sets";

operator record DataRecord "Coefficient data record for chemical definitions based on NASA source"
  extends Modelica.Icons.Record;

  Modelica.Units.SI.MolarMass MM "Molar mass";
  Modelica.Units.SI.MolarEnthalpy Hf "Enthalpy of formation at 298.15K, 1bar";
  Modelica.Units.SI.MolarEnthalpy H0 "H0(298.15K, 1bar) - H0(0K, 1bar)";
  Real alow[7] "Low temperature coefficients a at 298.15K, 1bar";
  Real blow[2] "Low temperature constants b at 298.15K, 1bar";
  Real ahigh[7] "High temperature coefficients a at 298.15K, 1bar";
  Real bhigh[2] "High temperature constants b at 298.15K, 1bar";

  Modelica.Units.SI.ChargeNumberOfIon z
  "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

  Chemical.Interfaces.Phase phase "State of matter";

  // following is used only if phase is not Chemical.Interfaces.Phase.Gas:
  Modelica.Units.SI.MolarVolume VmBase "Base molar volume at 298.15K, 1bar (molar volume Vm = VmBase + VmExcess, density = MM/Vm)";
  Modelica.Units.SI.MolarVolume VmExcess "Excess of molar volume at 298.15K, 1bar (activity coeficient = exp(VmExcess/VmBase)";

   encapsulated operator 'constructor'
    import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    import DataRecord=Chemical.Interfaces.DataRecord;
    import PhaseType=Chemical.Interfaces.Phase;

    constant Real R=1.380649e-23*6.02214076e23;
    constant Real T0=298.15 "Base temperature";
    constant Real p0=100000 "Base pressure";

     function fromModelicaDataRecord
      input ModelicaDataRecord mdata "Mass based data record";
      input Real z=0 "Charge number";
      input PhaseType phase=PhaseType.Gas "State of matter";
      input Real Vm=if (phase == PhaseType.Gas) then R*T0/p0 else 0.001*mdata.MM "Molar volume";
      input Real gamma=1 "Activity coefficient";

      output DataRecord result(
                  MM=mdata.MM,
                  Hf=mdata.Hf*mdata.MM,
                  H0=mdata.H0*mdata.MM,
                  alow=(mdata.R_s*mdata.MM/R) * mdata.alow,
                  blow=(mdata.R_s*mdata.MM/R) * mdata.blow,
                  ahigh=(mdata.R_s*mdata.MM/R) * mdata.ahigh,
                  bhigh=(mdata.R_s*mdata.MM/R) * mdata.bhigh,
                  z=z,
                  phase=phase,
                  VmBase=Vm/(1+log(gamma)),
                  VmExcess=Vm*log(gamma)/(1+log(gamma)))
                           "Molar based data record";
     algorithm
       assert(mdata.Tlimit==1000, "Tlimit must be 1000K!");
      annotation (Inline=true);
     end fromModelicaDataRecord;


    function fromValues
      input Real MM;
      input Real Hf;
      input Real H0;
      input Real alow[7];
      input Real blow[2];
      input Real ahigh[7];
      input Real bhigh[2];
      input Real z=0 "Charge number";
      input PhaseType phase=PhaseType.Incompressible "State of matter";
      input Real VmBase=if (phase == PhaseType.Gas) then R*T0/p0 else 0.001*MM "Base molar volume";
      input Real VmExcess=0 "Excess molar volume";

      output DataRecord result(
          MM = MM,
          Hf = Hf,
          H0 = H0,
          alow = alow,
          blow = blow,
          ahigh = ahigh,
          bhigh = bhigh,
          z = z,
          phase=phase,
          VmBase=VmBase,
          VmExcess=VmExcess);
    algorithm
      annotation (Inline=true);
    end fromValues;

   end 'constructor';
  annotation (Documentation(info="<html>
<p>
This data record contains the coefficients for the
ideal gas equations according to:
</p>
<blockquote>
  <p>McBride B.J., Zehe M.J., and Gordon S. (2002): <strong>NASA Glenn Coefficients
  for Calculating Thermodynamic Properties of Individual Species</strong>. NASA
  report TP-2002-211556</p>
</blockquote>
<p>
The equations have the following structure:
</p>
<div><img src=\"modelica://Modelica/Resources/Images/Media/IdealGases/Common/singleEquations.png\"></div>
<p>
The polynomials for h(T) and s0(T) are derived via integration from the one for cp(T)  and contain the integration constants b1, b2 that define the reference specific enthalpy and entropy. For entropy differences the reference pressure p0 is arbitrary, but not for absolute entropies. It is chosen as 1 standard atmosphere (101325 Pa).
</p>
<p>
For most gases, the region of validity is from 200 K to 6000 K.
The equations are split into two regions that are separated
by Tlimit (usually 1000 K). In both regions the gas is described
by the data above. The two branches are continuous and in most
gases also differentiable at Tlimit.
</p>
</html>"));
end DataRecord;

  connector ForeOld "Undirected connector outputting the forward state"

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

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.OutputSubstanceState state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.InputSubstanceState state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.OutputSolutionState solution "State of solution";
    stateOfMatter.OutputSubstanceData definition "Definition of substance";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={194,138,221},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid)}),
                              Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Fore is a output. </u>
</html>"));

  end ForeOld;

  connector RearOld "Undirected connector outputting the rearward state"

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

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    Chemical.Interfaces.InputSubstanceState state_forwards "State of substance in forwards direction";
    Chemical.Interfaces.OutputSubstanceState state_rearwards "State of substance in rearwards direction";

    Chemical.Interfaces.InputSolutionState solution "State of solution";
    stateOfMatter.InputSubstanceData definition "Definition of substance";

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={158,66,200},
            lineThickness=0.5,
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={194,138,221},
            lineThickness=0.5,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
                              Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<u>One of the two undirected connectors. The state information flows in both directions, forward and backward. </u>
<u>At positive molarflow, Rear is a input.</u>
</html>"));
  end RearOld;

 replaceable record SubstanceState "Set that defines a state of substance"
  extends Modelica.Icons.Record;

   Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";
   Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy of the substance";
 end SubstanceState;


 operator record SolutionState "Set that defines a state of solution"
  extends Modelica.Icons.Record;

   Modelica.Units.SI.Temperature T "Temperature of the solution";
   Modelica.Units.SI.Pressure p "Pressure of the solution";
   Modelica.Units.SI.ElectricPotential v "Electric potential in the solution";
   Modelica.Units.SI.AmountOfSubstance n "Amount of the solution";
   Modelica.Units.SI.Mass m "Mass of the solution";
   Modelica.Units.SI.Volume V "Volume of the solution";
   Modelica.Units.SI.Energy G "Free Gibbs energy of the solution";
   Modelica.Units.SI.ElectricCharge Q "Electric charge of the solution";
   Modelica.Units.SI.MoleFraction I "Mole fraction based ionic strength of the solution";

   encapsulated operator 'constructor'
     import Chemical.Interfaces.SolutionState;
     import Chemical.Interfaces.SolutionPort;
     import Chemical.Interfaces.Phase;
     constant Real R=1.380649e-23*6.02214076e23;
     constant Real T0=298.15 "Base temperature";
     constant Real p0=100000 "Base pressure";

     function fromValues
       input Phase phase "Phase of the chemical solution";
       input Real T=T0 "Temperature of the solution";
       input Real p=p0 "Pressure of the solution";
       input Real v=0 "Electric potential in the solution";
       input Real n=1 "Amount of the solution";
       input Real m=1 "Mass of the solution";
       input Real V=if (phase==Phase.Gas) then n*R*T/p else 1 "Volume of the solution";
       input Real G=0 "Free Gibbs energy of the solution";
       input Real Q=0 "Electric charge of the solution";
       input Real I=0 "Mole fraction based ionic strength of the solution";
       output SolutionState result(T=T,p=p,v=v,n=n,m=m,V=V,G=G,Q=Q,I=I);
     algorithm
       annotation(Inline = true);
     end fromValues;

     function fromSolutionPort
       input SolutionPort s;
       output SolutionState result(T=s.T,p=s.p,v=s.v,n=s.n,m=s.m,V=s.V,G=s.G,Q=s.Q,I=s.I);
     algorithm
       annotation(Inline = true);
     end fromSolutionPort;
   end 'constructor';

   encapsulated operator function '=='
     import Chemical.Interfaces.SolutionState;
     input SolutionState s1;
     input SolutionState s2;
     output Boolean result "= s1 == s2";
   algorithm
      result := s1.T == s2.T and s1.p == s2.p and s1.v == s2.v and s1.n == s2.n and s1.m == s2.m and s1.V == s2.V and s1.G == s2.G and s1.Q == s2.Q and s1.I == s2.I;
      annotation(Inline = true);
   end '==';
 end SolutionState;

 operator record SolutionStateParameters "Set that defines a state of solution"
  extends Modelica.Icons.Record;

   parameter Modelica.Units.SI.Temperature T=293.15   "Temperature of the solution";
   parameter Modelica.Units.SI.Pressure p=101325   "Pressure of the solution";
   parameter Modelica.Units.SI.Mass m = 1 "Mass of the solution";
   parameter Modelica.Units.SI.Volume V(displayUnit="l")=0.001  "Volume of the solution";
   parameter Modelica.Units.SI.AmountOfSubstance n = 1 "Amount of the solution";
   parameter Modelica.Units.SI.ElectricPotential v = 0 "Electric potential in the solution";
   parameter Modelica.Units.SI.Energy G = 0 "Free Gibbs energy of the solution";
   parameter Modelica.Units.SI.ElectricCharge Q = 0 "Electric charge of the solution";
   parameter Modelica.Units.SI.MoleFraction I = 0 "Mole fraction based ionic strength of the solution";


 end SolutionStateParameters;


  connector InputDefinition = input Chemical.Interfaces.Definition;
  connector InputSubstanceDefinition
                            = input Chemical.Interfaces.SubstanceDefinition;
  connector InputSubstanceState = input SubstanceState;
  connector InputSolutionState = input SolutionState;
  connector OutputDefinition = output Chemical.Interfaces.Definition;
  connector OutputSubstanceDefinition
                             = output Chemical.Interfaces.Definition;
  connector OutputSubstanceState = output SubstanceState;
  connector OutputSolutionState = output SolutionState;

  partial package StateOfMatter "Abstract package for all state of matters"



   replaceable partial record SubstanceDataParameters
      "Definition data of the chemical substance as parameters"

   end SubstanceDataParameters;

   replaceable record SubstanceData "Minimal set that defines a substance"
    extends Modelica.Icons.Record;


   end SubstanceData;

    replaceable connector InputSubstanceData

    end InputSubstanceData;

    replaceable connector OutputSubstanceData

    end OutputSubstanceData;



    replaceable partial model BaseProperties "Base properties of the substance"

      parameter Boolean FixedSubstanceData "substanceDataVar==substanceData";

      parameter SubstanceDataParameters substanceData "used only if FixedSubstanceData to help initialization";

      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

      parameter Boolean SolutionObserverOnly = false "True if substance does not affect the solution";


      InputSubstanceData substanceDataVar "Definition of the substance";

      InputSolutionState solutionState "State of the solution";

      InputAmountOfSubstance amountOfBaseMolecules
        "Amount of base molecules inside all clusters in compartment";

      InputMolarFlowRate n_flow "Molar change of the substance";

      InputHeatFlowRate h_flow "Substance enthaply change";


      Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
        "Molar concentration of particles";

      Modelica.Units.SI.MassConcentration M(displayUnit="mg/l")
            "Mass concentration";

      Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
            "Molality";

      Modelica.Units.SI.MassFraction X "Mass fraction";


      Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";

      Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";


      Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of the substance";

      Modelica.Units.SI.MolarEnthalpy hPure
      "Molar enthalpy of the pure substance";

      Modelica.Units.SI.MolarEntropy sPure
      "Molar entropy of the pure substance";

      Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

      Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

      Modelica.Units.SI.MolarVolume Vm
      "Molar volume of the substance";

      Modelica.Units.SI.MolarMass MM
      "Molar mass of the substance";

      Modelica.Units.SI.MolarVolume VmPure
      "Molar volume of the pure substance";

      Modelica.Units.SI.MolarVolume VmExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

     // Local connector definition, used for equation balancing check

      connector InputMolarFlowRate = input Modelica.Units.SI.MolarFlowRate
        "Molar flow rate as input signal connector";
      connector InputHeatFlowRate = input Modelica.Units.SI.HeatFlowRate
        "Heat flow rate as input signal connector";

      connector InputSubstanceData = input SubstanceData
        "Substance definition as input signal connector";
      connector InputTemperature = input Modelica.Units.SI.Temperature
        "Temperature as input signal connector";
      connector InputPressure = input Modelica.Units.SI.Pressure
        "Pressure as input signal connector";
      connector InputElectricPotential = input
          Modelica.Units.SI.ElectricPotential
        "Electric potential as input signal connector";
      connector InputMoleFraction = input Modelica.Units.SI.MoleFraction
        "Mole fraction as input signal connector";
      connector InputAmountOfSubstance = input
          Modelica.Units.SI.AmountOfSubstance
        "Amount of substance as input signal connector";





      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*specificAmountOfFreeBaseMolecule(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*specificAmountOfParticles(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solutionState.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          selfClusteringBondEnthalpy(substanceDataVar)
        - solutionState.T * selfClusteringBondEntropy(substanceDataVar)
        "Gibbs energy of hydrogen bond between H2O molecules";


      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";


      outer Modelica.Fluid.System system "System wide properties";

      Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;

    equation
      assert(x > 0, "Molar fraction must be positive");

     gamma = activityCoefficient(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     z = chargeNumberOfIon(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);

     MM = 1/specificAmountOfParticles(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     h = molarEnthalpy(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     sPure = molarEntropyPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     hPure = uPure + solutionState.T*sPure;
     u0 = chemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     uPure = electroChemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     Vm = molarVolume(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     VmPure = molarVolumePure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     VmExcess = molarVolumeExcess(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     //  molarHeatCapacityCp = smolarHeatCapacityCp(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     u = chemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I)
       + (Modelica.Constants.R*solutionState.T)*log(a)
       + z*Modelica.Constants.F*solutionState.v;



      // during initialization it does not take value from parameter substanceData if not useInlet,
      // so instead of just "selfClustering(substanceDataVar)" it must be written
      if (not FixedSubstanceData and selfClustering(substanceData)) or selfClustering(substanceDataVar) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solutionState.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solutionState.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solutionState.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        x = ((2*SelfClustering_K+solutionState.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solutionState.n/amountOfBaseMolecules)+(solutionState.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties

       //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
       //enthalpy = h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        dH = der(h)*
          amountOfBaseMolecules + n_flow*h - h_flow +
          selfClusteringBondEnthalpy(substanceDataVar)*der(amountOfBonds)
                        "heat transfer from other substances in solution [J/s]";

        Gj =amountOfBaseMolecules*u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
        //enthalpy = h*amountOfBaseMolecules;
        dH =  der(h)*amountOfBaseMolecules + n_flow*h - h_flow
                  "heat transfer from other substances in solution [J/s]";

        Gj = amountOfBaseMolecules*u "Gibbs energy of the substance [J]";

      end if;


      x = amountOfFreeMolecule/solutionState.n "mole fraction [mol/mol]";

      c = amountOfParticles/solutionState.V "concentration [mol/m3]";

      M = amountOfParticles*molarMassOfBaseMolecule(substanceDataVar)/solutionState.V "mass concentration [kg/m3]";

      b = amountOfParticles/solutionState.m "molality [mol/kg]";

      X  = amountOfParticles*molarMassOfBaseMolecule(substanceDataVar)/solutionState.m "mass fraction [kg/kg]";

      if SolutionObserverOnly then
        i = 0;
        dV = 0;
        nj = 0;
        mj = 0;
        Vj = 0;
        Qj = 0;
        Ij = 0;
      else
      //solution flows
        i = Modelica.Constants.F*z*n_flow +
          Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
        dV = Vm*n_flow + der(Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
        nj = amountOfParticles;
        mj = amountOfBaseMolecules*molarMassOfBaseMolecule(substanceDataVar);
        Vj = amountOfBaseMolecules*Vm;
        Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
        Ij = (1/2)*(amountOfBaseMolecules*z^2);
      end if;


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end BaseProperties;



   replaceable function processData
     "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
        extends Modelica.Icons.Function;
    input Real K "Process dissociation constant (mole-fraction based) at 25°C,1bar";
    input Modelica.Units.SI.MolarEnergy dH "Process molar enthalpy change at 25°C,1bar";
    input Modelica.Units.SI.MolarHeatCapacity dCp "Process molar heat capacity change at 25°C,1bar";
    input Modelica.Units.SI.SpecificVolume dVs "Process specific volume change at 25°C,1bar";
    output SubstanceData processData "Data record of process changes";
   end processData;

   replaceable function firstProductDefinition
       "Return formation definition of the first product of chemical process"
        extends Modelica.Icons.Function;
    //input
    input Modelica.Units.SI.StoichiometricNumber s[:] "Stoichiometric reaction coefficient for substrates [nS]";
    input Modelica.Units.SI.StoichiometricNumber p[:] "Stoichiometric reaction coefficient for products [nP]";
    input SubstanceData processData "Data record of process changes";
    input SubstanceData substratesData[:] "Substrates definitions [nS]";
    input SubstanceData productsData[:] "Other products definitions [nP-1]";
    output SubstanceData firstProductDefinition "Definition of the first product in process";
   end firstProductDefinition;

  /*
 replaceable function substrateFromIncompressible
    "Cast substrate definition from Incompressible (but the phase remains Incompressinble) - use only for firstProductDefinition, where phase shift will be included"
     extends Modelica.Icons.Function;
     input Chemical.Interfaces.Incompressible.SubstanceData dataIncompressible "Substance definition as Incompressible";
     output SubstanceData data "Casted substrate definition";
 end substrateFromIncompressible;

  replaceable function substrateFromIdealGas
    "Cast substrate definition from IdealGas (but the phase remains IdealGas) - use only for firstProductDefinition, where phase shift will be included"
     extends Modelica.Icons.Function;
     input Chemical.Interfaces.IdealGas.SubstanceData dataIdealGas "Substance definition as IdealGas";
     output SubstanceData data "Casted substrate definition";
  end substrateFromIdealGas;

  replaceable function substrateFromIdealGasMSL
    "Cast substrate definition from IdealGas (but the phase remains IdealGasMSL) - use only for firstProductDefinition, where phase shift will be included"
     extends Modelica.Icons.Function;
     input Chemical.Interfaces.IdealGasMSL.SubstanceData dataIdealGasMSL "Substance definition as IdealGasMSL";
     output SubstanceData data "Casted substrate definition";
  end substrateFromIdealGasMSL;

  replaceable function substrateFromIdealGasShomate
    "Cast substrate definition from IdealGas (but the phase remains IdealGasShomate) - use only for firstProductDefinition, where phase shift will be included"
     extends Modelica.Icons.Function;
     input Chemical.Interfaces.IdealGasShomate.SubstanceData dataIdealGasShomate "Substance definition as IdealGasShomate";
     output SubstanceData data "Casted substrate definition";
  end substrateFromIdealGasShomate;
*/

   replaceable function activityCoefficient
    "Return activity coefficient of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

      output Real activityCoefficient "Activity Coefficient";
   end activityCoefficient;

   replaceable function chargeNumberOfIon
    "Return charge number of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.ChargeNumberOfIon chargeNumberOfIon
      "Charge number of ion";
   end chargeNumberOfIon;

   replaceable function molarEnthalpyElectroneutral
    "Molar enthalpy of the substance in electroneutral solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEnthalpy molarEnthalpyElectroneutral
      "Molar enthalpy";
   end molarEnthalpyElectroneutral;

   function molarEnthalpy
    "Molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEnthalpy molarEnthalpy
      "Molar enthalpy";
   algorithm
      molarEnthalpy := molarEnthalpyElectroneutral(substanceData,T,p,v,I) +
           Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
      annotation (Inline=true, smoothOrder=2);
   end molarEnthalpy;

   replaceable function molarEntropyPure
    "Molar entropy of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEntropy molarEntropyPure
      "Molar entropy of the pure substance";
   end molarEntropyPure;

    function molarEntropy "Molar entropy of the substance in the solution"
          extends Modelica.Icons.Function;
    input Modelica.Units.SI.ChemicalPotential u
      "Electro-chemical potential of the substance";
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarEntropy molarEntropy "Molar entropy";
    algorithm
        molarEntropy :=  (u - molarEnthalpy(substanceData,T,p,v,I))/T;
    end molarEntropy;

   function chemicalPotentialPure "Chemical potential of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.ChemicalPotential chemicalPotentialPure
      "Base chemical potential";
   algorithm
       chemicalPotentialPure :=  molarEnthalpyElectroneutral(substanceData,T,p,v,I) - T*molarEntropyPure(substanceData,T,p,v,I);
   end chemicalPotentialPure;

   function electroChemicalPotentialPure
    "Electro-chemical potential of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.ChemicalPotential
      electroChemicalPotentialPure "Base electro-chemical potential";
   algorithm
    electroChemicalPotentialPure := chemicalPotentialPure(
         substanceData,
         T,
         p,
         v,
         I) + Modelica.Constants.F*chargeNumberOfIon(substanceData,T,p,v,I)*v;
   end electroChemicalPotentialPure;

   replaceable function molarVolumePure "Molar volume of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
   end molarVolumePure;

   function molarVolumeExcess
    "Excess molar volume of the substance in the solution"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarVolume molarVolumeExcess
      "Excess molar volume of the substance in the solution";
   algorithm
      molarVolumeExcess := molarVolumePure(substanceData,T,p,v,I)*
         log(activityCoefficient(substanceData,T,p,v,I)); //zero if activityCoefficient==1
      annotation (Inline=true, smoothOrder=2);
   end molarVolumeExcess;

   replaceable function molarVolume "Molar volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.MolarVolume molarVolume "Molar volume";
   algorithm
    molarVolume :=molarVolumePure(
         substanceData,
         T,
         p,
         v,
         I) + molarVolumeExcess(
         substanceData,
         T,
         p,
         v,
         I);
      annotation (Inline=true, smoothOrder=2);
   end molarVolume;

   replaceable function molarHeatCapacityCp
    "Molar heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.MolarHeatCapacity molarHeatCapacityCp
      "Molar heat capacity at constant pressure";
   end molarHeatCapacityCp;

   replaceable function molarMassOfBaseMolecule
      "Molar mass of base molecule of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarMass molarMass "Molar mass";
   end molarMassOfBaseMolecule;

   replaceable function selfClustering "returns true if substance molecules are joining together to clusters"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
          output Boolean selfClustering;
   algorithm
     selfClustering:=false;
   end selfClustering;

   replaceable function selfClusteringBondEnthalpy
    "Enthalpy of joining two base molecules of the substance together to cluster"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
   algorithm
     selfClusteringEnthalpy:=0;
   end selfClusteringBondEnthalpy;

   replaceable function selfClusteringBondEntropy
    "Entropy of joining two base molecules of the substance together to cluster"
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
   algorithm
     selfClusteringEntropy:=0;
   end selfClusteringBondEntropy;

   replaceable function selfClusteringBondVolume
       extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarVolume selfClusteringBondVolume;
   algorithm
     selfClusteringBondVolume:=0;
   end selfClusteringBondVolume;

   replaceable function selfClusteringBondHeatCapacityCp
      extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    output Modelica.Units.SI.MolarHeatCapacity selfClusteringBondHeatCapacityCp;
   algorithm
     selfClusteringBondHeatCapacityCp:=0;
   end selfClusteringBondHeatCapacityCp;

    replaceable function specificAmountOfParticles
      "Amount of particles per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg")
        "Amount of substance particles per its mass";
    algorithm
      specificAmountOfSubstance := 1/molarMassOfBaseMolecule(substanceData);
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfParticles;

    replaceable function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";
    algorithm
      specificAmountOfFreeBaseMolecule := 1/molarMassOfBaseMolecule(substanceData);

      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

  /* replaceable function solution_temperature_
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
     extends Modelica.Icons.Function;
    input SubstanceData substanceData[:] "Data record of substances";
  input Modelica.Units.SI.MolarEnthalpy h
    "Molar enthalpy of solution (x*substances_h)";
  input Modelica.Units.SI.MoleFraction x[:]
    "Mole fractions of substances";
  input Modelica.Units.SI.Pressure p=100000 "Pressure";
  input Modelica.Units.SI.ElectricPotential v=0
    "Electric potential of the substance";
  input Modelica.Units.SI.MoleFraction I=0
    "Ionic strengh (mole fraction based)";

  output Modelica.Units.SI.Temperature T "Temperature";
    annotation (__Dymola_DymolaStoredErrors(thetext="/*replaceable function solution_temperature_
  \"Temperature of the solution from specific enthalpy and mass fractions of substances\"
    extends Modelica.Icons.Function;
   input SubstanceData substanceData[:] \"Data record of substances\";
 input Modelica.Units.SI.MolarEnthalpy h
   \"Molar enthalpy of solution (x*substances_h)\";
 input Modelica.Units.SI.MoleFraction x[:]
   \"Mole fractions of substances\";
 input Modelica.Units.SI.Pressure p=100000 \"Pressure\";
 input Modelica.Units.SI.ElectricPotential v=0
   \"Electric potential of the substance\";
 input Modelica.Units.SI.MoleFraction I=0
   \"Ionic strengh (mole fraction based)\";

 output Modelica.Units.SI.Temperature T \"Temperature\";
"));
end solution_temperature_;
*/

   replaceable function specificEnthalpy
     "Specific molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
      "Specific enthalpy";

   algorithm

     specificEnthalpy := molarEnthalpy(
       substanceData,
       T,
       p,
       v,
       I)/
       molarMassOfBaseMolecule(substanceData);
   end specificEnthalpy;

   replaceable function specificVolume "Specific volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.SpecificVolume specificVolume "Specific volume";

   algorithm

    specificVolume := molarVolume(
         substanceData,
         T,
         p,
         v,
         I) /
       molarMassOfBaseMolecule(substanceData);
   end specificVolume;

    replaceable function specificHeatCapacityCp
    "Specific heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";
    output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
      "Specific heat capacity at constant pressure";

    algorithm

    specificHeatCapacityCp := molarHeatCapacityCp(
         substanceData,
         T,
         p,
         v,
         I) /
       molarMassOfBaseMolecule(substanceData);
    end specificHeatCapacityCp;

   replaceable function temperature
    "Temperature of the substance from its enthalpy"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";
   end temperature;

   replaceable function solution_temperature
    "Temperature of the solution from specific enthalpy and mass fractions of substances"
       extends Modelica.Icons.Function;
      input SubstanceData substanceData[:] "Data record of substances";
    input Modelica.Units.SI.SpecificEnthalpy h
      "Specific enthalpy of solution";
    input Modelica.Units.SI.MassFraction X[:]
      "Mass fractions of substances";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Temperature T "Temperature";
   end solution_temperature;

   replaceable function density
        "Return density of the substance in the solution"
          extends Modelica.Icons.Function;
          input SubstanceData substanceData "Data record of substance";
    input Modelica.Units.SI.Temperature T=298.15 "Temperature";
    input Modelica.Units.SI.Pressure p=100000 "Pressure";
    input Modelica.Units.SI.ElectricPotential v=0
      "Electric potential of the substance";
    input Modelica.Units.SI.MoleFraction I=0
      "Ionic strengh (mole fraction based)";

    output Modelica.Units.SI.Density density "Density";
   end density;

    replaceable partial model BaseProperties2 "Base properties of the substance"

      parameter Boolean FixedSubstanceData "substanceDataVar==substanceData";

      parameter SubstanceDataParameters substanceData "used only if FixedSubstanceData to help initialization";

      parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

      parameter Boolean SolutionObserverOnly = false "True if substance does not affect the solution";

      InputSubstanceData substanceDataVar "Definition of the substance";

      InputSolutionState solutionState "State of the solution";

      InputAmountOfSubstance amountOfBaseMolecules
        "Amount of base molecules inside all clusters in compartment";

      InputMolarFlowRate n_flow "Molar change of the substance";

      InputHeatFlowRate h_flow "Substance enthaply change";

      Modelica.Units.SI.MoleFraction x "Mole fraction of the substance";

      Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
        "Molar concentration of particles";

      Modelica.Units.SI.MassConcentration M(displayUnit="mg/l")
            "Mass concentration";

      Modelica.Units.SI.Molality b(displayUnit="mmol/kg")
            "Molality";

      Modelica.Units.SI.MassFraction X "Mass fraction";

      Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential of the substance";

      Modelica.Units.SI.ActivityOfSolute a
      "Activity of the substance (mole-fraction based)";

      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of ion";

      Modelica.Units.SI.MolarEnthalpy h
      "Molar enthalpy of the substance";

      Modelica.Units.SI.MolarEnthalpy hPure
      "Molar enthalpy of the pure substance";

      Modelica.Units.SI.MolarEntropy sPure
      "Molar entropy of the pure substance";

      Modelica.Units.SI.ChemicalPotential u0
      "Chemical potential of the pure substance";

      Modelica.Units.SI.ChemicalPotential uPure
      "Electro-Chemical potential of the pure substance";

      Modelica.Units.SI.MolarVolume Vm
      "Molar volume of the substance";

      Modelica.Units.SI.MolarMass MM
      "Molar mass of the substance";

      Modelica.Units.SI.MolarVolume VmPure
      "Molar volume of the pure substance";

      Modelica.Units.SI.MolarVolume VmExcess
      "Molar volume excess of the substance in solution (typically it is negative as can be negative)";

     // Local connector definition, used for equation balancing check

      connector InputMolarFlowRate = input Modelica.Units.SI.MolarFlowRate
        "Molar flow rate as input signal connector";
      connector InputHeatFlowRate = input Modelica.Units.SI.HeatFlowRate
        "Heat flow rate as input signal connector";

      connector InputSubstanceData = input SubstanceData
        "Substance definition as input signal connector";
      connector InputTemperature = input Modelica.Units.SI.Temperature
        "Temperature as input signal connector";
      connector InputPressure = input Modelica.Units.SI.Pressure
        "Pressure as input signal connector";
      connector InputElectricPotential = input
          Modelica.Units.SI.ElectricPotential
        "Electric potential as input signal connector";
      connector InputMoleFraction = input Modelica.Units.SI.MoleFraction
        "Mole fraction as input signal connector";
      connector InputAmountOfSubstance = input
          Modelica.Units.SI.AmountOfSubstance
        "Amount of substance as input signal connector";

      Modelica.Units.SI.AmountOfSubstance amountOfFreeMolecule(start=
           m_start*specificAmountOfFreeBaseMolecule(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of free molecules not included inside any clusters in compartment";

      Modelica.Units.SI.AmountOfSubstance amountOfParticles(start=
           m_start*specificAmountOfParticles(
                                       substanceData,
                                       T=system.T_ambient,
                                       p=system.p_ambient))
        "Amount of particles/clusters in compartment";

      Modelica.Units.SI.MoleFraction SelfClustering_K=exp(-SelfClustering_dG/(
          Modelica.Constants.R*solutionState.T))
        "Dissociation constant of hydrogen bond between base molecules";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dG=
          selfClusteringBondEnthalpy(substanceDataVar)
        - solutionState.T * selfClusteringBondEntropy(substanceDataVar)
        "Gibbs energy of hydrogen bond between H2O molecules";

      Modelica.Units.SI.AmountOfSubstance amountOfBonds
        "Amount of hydrogen bonds between molecules in compartment";

      outer Modelica.Fluid.System system "System wide properties";

      Real i,dH,dV,nj,mj,Vj,Gj,Qj,Ij;

    equation
      assert(x > 0, "Molar fraction must be positive");

     gamma = activityCoefficient(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     z = chargeNumberOfIon(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);

     MM = 1/specificAmountOfParticles(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     h = molarEnthalpy(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     sPure = molarEntropyPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     hPure = uPure + solutionState.T*sPure;
     u0 = chemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     uPure = electroChemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     Vm = molarVolume(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     VmPure = molarVolumePure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     VmExcess = molarVolumeExcess(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);
     //  molarHeatCapacityCp = smolarHeatCapacityCp(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I);

     //activity of the substance
     a = gamma*x;

     //electro-chemical potential of the substance in the solution
     u = chemicalPotentialPure(substanceDataVar,solutionState.T,solutionState.p,solutionState.v,solutionState.I)
       + (Modelica.Constants.R*solutionState.T)*log(a)
       + z*Modelica.Constants.F*solutionState.v;

      // during initialization it does not take value from parameter substanceData if not useInlet,
      // so instead of just "selfClustering(substanceDataVar)" it must be written
      if (not FixedSubstanceData and selfClustering(substanceData)) or selfClustering(substanceDataVar) then

        //Liquid cluster theory - equilibrium:
        //x[i] = x*(K*x)^i .. mole fraction of cluster composed with i base molecules
        //amountOfParticles/solutionState.n = x/(1-K*x);                //sum(x[i])
        //amountOfBaseMolecules/solutionState.n = x/((1-K*x)^2);            //sum(i*x[i])
        //amountOfHydrogenBonds/solutionState.n = x*x*K/((1-K*x)^2);   //sum((i-1)*x[i])

        amountOfParticles*(1 - SelfClustering_K*x) = amountOfFreeMolecule;

        //Calculation of "abs(amountOfBaseMolecules*(1 - SelfClustering_K*x)) = amountOfParticles":
        x = ((2*SelfClustering_K+solutionState.n/amountOfBaseMolecules) - sqrt((4*SelfClustering_K*solutionState.n/amountOfBaseMolecules)+(solutionState.n/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfBonds = amountOfBaseMolecules*x*SelfClustering_K;

        //TODO: may be the volume of the same number of free water molecules is different as volume of the same number of water molecules in cluster ..
        //TODO: more precise calculation of other properties

       //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
       //enthalpy = h*amountOfBaseMolecules + amountOfAdditionalBonds*bondEnthalpy;
        dH = der(h)*
          amountOfBaseMolecules + n_flow*h - h_flow +
          selfClusteringBondEnthalpy(substanceDataVar)*der(amountOfBonds)
                        "heat transfer from other substances in solution [J/s]";

        Gj =amountOfBaseMolecules*u + amountOfBonds*SelfClustering_dG
                        "Gibbs energy of the substance";

      else

        amountOfParticles = amountOfFreeMolecule;
        amountOfBaseMolecules = amountOfFreeMolecule;
        amountOfBonds = 0;

        //der(enthalpy) = dH + n_flow*actualStream(port_a.h_outflow);
        //enthalpy = h*amountOfBaseMolecules;
        dH =  der(h)*amountOfBaseMolecules + n_flow*h - h_flow
                  "heat transfer from other substances in solution [J/s]";

        Gj = amountOfBaseMolecules*u "Gibbs energy of the substance [J]";

      end if;

      x = amountOfFreeMolecule/solutionState.n "mole fraction [mol/mol]";

      c = amountOfParticles/solutionState.V "concentration [mol/m3]";

      M = amountOfParticles*molarMassOfBaseMolecule(substanceDataVar)/solutionState.V "mass concentration [kg/m3]";

      b = amountOfParticles/solutionState.m "molality [mol/kg]";

      X  = amountOfParticles*molarMassOfBaseMolecule(substanceDataVar)/solutionState.m "mass fraction [kg/kg]";

      if SolutionObserverOnly then
        i = 0;
        dV = 0;
        nj = 0;
        mj = 0;
        Vj = 0;
        Qj = 0;
        Ij = 0;
      else
      //solution flows
        i = Modelica.Constants.F*z*n_flow +
          Modelica.Constants.F*der(z)*amountOfBaseMolecules "change of sunstance charge [A]";
        dV = Vm*n_flow + der(Vm)*amountOfBaseMolecules "change of substance volume [m3/s]";

      //extensive properties
        nj = amountOfParticles;
        mj = amountOfBaseMolecules*molarMassOfBaseMolecule(substanceDataVar);
        Vj = amountOfBaseMolecules*Vm;
        Qj = Modelica.Constants.F*amountOfBaseMolecules*z;
        Ij = (1/2)*(amountOfBaseMolecules*z^2);
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end BaseProperties2;
    annotation (Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end StateOfMatter;

  package Incompressible "Incompressible as basic state of matter"
    extends StateOfMatter;

    constant String StateName="Incompressible"; //type checking
    constant Modelica.Units.SI.Temperature T0=298.15;
    constant Modelica.Units.SI.Pressure p0=100000;

    redeclare record extends SubstanceDataParameters "Base substance data"

      parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=1 "Molar weight of the substance";

      parameter Modelica.Units.SI.ChargeNumberOfIon z=0
        "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=0
        "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=0
        "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.ActivityCoefficient gamma=1
        "Activity coefficient of the substance";

      parameter Modelica.Units.SI.MolarHeatCapacity Cp=1
        "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";

      parameter Real SelfClustering=0
        "Pure substance is making clusters (weak bonds between molecules)";
                                      //false

      parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
        "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";
      //-20000
      parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
        "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

      parameter Modelica.Units.SI.SpecificVolume Vs(displayUnit="dm3/kg") = 0.001
        "Specific volume of the pure substance (default 1L/kg)";

      //      parameter Modelica.SIunits.MolarHeatCapacity Cv = Cp
      //      "Molar heat capacity of the substance at constant volume";

      annotation (preferredView="info", Documentation(revisions="<html>
<p><i>2015-2025</i></p>
<p>Marek Mateják </p>
</html>"));
    end SubstanceDataParameters;

    redeclare operator record extends SubstanceData "Base substance data"

      Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa") "Molar weight of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol") "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol") "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      Modelica.Units.SI.ActivityCoefficient gamma "Activity coefficient of the substance";

      Modelica.Units.SI.MolarHeatCapacity Cp "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";

      Real SelfClustering "Pure substance is making clusters (weak bonds between molecules)";

      Modelica.Units.SI.ChemicalPotential SelfClustering_dH "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";

      Modelica.Units.SI.MolarEntropy SelfClustering_dS "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";

      Modelica.Units.SI.SpecificVolume Vs(displayUnit="dm3/kg") "Specific volume of the pure substance (default 1L/kg)";

      // Modelica.Units.SI.Density density(displayUnit="kg/dm3") "Density of the pure substance (default density of water at 25degC)";

      encapsulated operator 'constructor'
        import Modelica.Units.SI.MolarMass;
        import SD = Chemical.Interfaces.Incompressible.SubstanceData;
        import IdealGasSD = Chemical.Interfaces.IdealGas.SubstanceData;
        import IdealGasMSLSD = Chemical.Interfaces.IdealGasMSL.SubstanceData;
        import IdealGasShomateSD = Chemical.Interfaces.IdealGasShomate.SubstanceData;
        import ReferenceEnthalpy = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;
        import Modelica.Media.IdealGases.Common.Functions.h_T;
        import Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral;
        import Chemical.Interfaces.IdealGasMSL.molarEntropyPure;
        import Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp;
        import Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral;
        import Chemical.Interfaces.IdealGasShomate.molarEntropyPure;
        import Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp;

        function fromReal
          input Real MolarWeight=1 "Molar weight of the substance";
          input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
          input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";
          input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";
          input Real gamma=1 "Activity coefficient of the substance";
          input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
          input Real SelfClustering=0 "Pure substance is making clusters (weak bonds between molecules)";
          input Real SelfClustering_dH=0 "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";
          input Real SelfClustering_dS=0 "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";
          input Real Vs=0.001 "Specific volume of the pure substance (default 1L/kg)";
          output SD data=SD(
                    MolarWeight=MolarWeight,
                    z=z,
                    DfG=DfG,
                    DfH=DfH,
                    gamma=gamma,
                    Cp=Cp,
                    SelfClustering=SelfClustering,
                    SelfClustering_dH=SelfClustering_dH,
                    SelfClustering_dS=SelfClustering_dS,
                    Vs=Vs);
        algorithm
          annotation (Inline=true);
        end fromReal;

        function fromIncompressible
          input SD dataIncompressible;
          output SD data=dataIncompressible;
        algorithm
          annotation (Inline=true);
        end fromIncompressible;

        function fromIdealGas
          input IdealGasSD dataIdealGas;
          output SD data=SD(
                    MolarWeight=dataIdealGas.MolarWeight,
                    z=dataIdealGas.z,
                    DfG=dataIdealGas.DfG,
                    DfH=dataIdealGas.DfH,
                    gamma=1,
                    Cp=dataIdealGas.Cp,
                    SelfClustering=0,
                    SelfClustering_dH=0,
                    SelfClustering_dS=0,
                    Vs=1);
          //TODO check specific volume - Vs
        algorithm
          annotation (Inline=true);

        end fromIdealGas;


        function fromIdealGasMSL
              input IdealGasMSLSD dataIdealGasMSL;
              output SD data=SD(
                        MolarWeight=dataIdealGasMSL.data.MM,
                        z=dataIdealGasMSL.z,
                        DfG=Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(
                          dataIdealGasMSL,
                          298.15,
                          100000,
                          0,
                          0) - 298.15*Chemical.Interfaces.IdealGasMSL.molarEntropyPure(
                          dataIdealGasMSL,
                          298.15,
                          100000,
                          0,
                          0),
                        DfH=Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(
                          dataIdealGasMSL,
                          298.15,
                          100000,
                          0,
                          0),
                        gamma=1,
                        Cp=Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp(
                          dataIdealGasMSL,
                          298.15,
                          100000,
                          0,
                          0),
                        SelfClustering=0,
                        SelfClustering_dH=0,
                        SelfClustering_dS=0,
                        Vs=1);
              //TODO check specific volume - Vs
        algorithm
              annotation (Inline=true);
        end fromIdealGasMSL;


        /*
  function fromIdealGasMSL
    input IdealGasMSLSD dataIdealGasMSL;
    output SD data;
    //TODO check specific volume - Vs
  protected 
    Real DfH,DfS,Cp;
  algorithm 

    DfH := Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(
          dataIdealGasMSL,
          298.15,
          100000,
          0,
          0);

    DfS := Chemical.Interfaces.IdealGasMSL.molarEntropyPure(
          dataIdealGasMSL,
          298.15,
          100000,
          0,
          0);

    Cp := Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp(
          dataIdealGasMSL,
          298.15,
          100000,
          0,
          0);

    data := SD(
              MolarWeight=dataIdealGasMSL.data.MM,
              z=dataIdealGasMSL.z,
              DfG=DfH - 298.15*DfS,
              DfH=DfH,
              gamma=1,
              Cp=Cp,
              SelfClustering=0,
              SelfClustering_dH=0,
              SelfClustering_dS=0,
              Vs=1);
   // annotation (Inline=true);
  end fromIdealGasMSL;

*/
        function fromIdealGasShomate
          input IdealGasShomateSD dataIdealGasShomate;
          output SD data=SD(
                    MolarWeight=dataIdealGasShomate.MolarWeight,
                    z=dataIdealGasShomate.z,
                    DfG=Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
                      dataIdealGasShomate,
                      298.15,
                      100000,
                      0,
                      0) - 298.15*Chemical.Interfaces.IdealGasShomate.molarEntropyPure(
                      dataIdealGasShomate,
                      298.15,
                      100000,
                      0,
                      0),
                    DfH=Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
                      dataIdealGasShomate,
                      298.15,
                      100000,
                      0,
                      0),
                    gamma=1,
                    Cp=Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(
                      dataIdealGasShomate,
                      298.15,
                      100000,
                      0,
                      0),
                    SelfClustering=0,
                    SelfClustering_dH=0,
                    SelfClustering_dS=0,
                    Vs=1);
          //TODO check specific volume - Vs
        algorithm
          annotation (Inline=true);
        end fromIdealGasShomate;


      end 'constructor';
      annotation (preferredView="info", Documentation(revisions="<html>
<p><i>2015-2025</i></p>
<p>Marek Mateják </p>
</html>"));
    end SubstanceData;


    redeclare connector InputSubstanceData = input SubstanceData;
    redeclare connector OutputSubstanceData = output SubstanceData;


    redeclare replaceable model extends BaseProperties "Base properties of incompressible substance"
    end BaseProperties;

    redeclare function extends processData
     "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
    algorithm
      processData :=  SubstanceDataParameters(
        MolarWeight = 0,
        z = 0,
        DfG = -Modelica.Constants.R*T0*log(K),
        DfH = dH,
        gamma = 1,
        Cp = dCp,
        SelfClustering = 0,
        SelfClustering_dH = 0,
        SelfClustering_dS = 0,
        Vs = dVs);
    end processData;

    redeclare function extends firstProductDefinition
       "Return formation definition of the first product of chemical process"
    algorithm
       firstProductDefinition :=  SubstanceDataParameters(
         MolarWeight = ((p[2:end]*productsData.MolarWeight) - (s*substratesData.MolarWeight))/p[1],
         z = ((p[2:end]*productsData.z) - (s*substratesData.z))/p[1],
         DfG = ((p[2:end]*productsData.DfG) - (s*substratesData.DfG) - processData.DfG)/p[1],
         DfH = ((p[2:end]*productsData.DfH) - (s*substratesData.DfH) - processData.DfH)/p[1],
         gamma = 1,
         Cp = ((p[2:end]*productsData.Cp) - (s*substratesData.Cp) - processData.Cp)/p[1],
         SelfClustering = 0,
         SelfClustering_dH = 0,
         SelfClustering_dS = 0,
         Vs = ((p[2:end]*productsData.Vs) - (s*substratesData.Vs) - processData.Vs)/p[1]);
    end firstProductDefinition;




   //  redeclare function extends substrateFromIncompressible
   //   "Cast substrate definition from Incompressible (but the phase remains Incompressinble) - use only for firstProductDefinition, where phase shift will be included"
   //  algorithm
   //    data := dataIncompressible;
        /*SubstanceData(
       MolarWeight = dataIncompressible.MolarWeight,
       z = dataIncompressible.z,
       DfG = dataIncompressible.DfG,
       DfH = dataIncompressible.DfH,
       gamma = dataIncompressible.gamma,
       Cp = dataIncompressible.Cp,
       SelfClustering = dataIncompressible.SelfClustering,
       SelfClustering_dH = dataIncompressible.SelfClustering_dH,
       SelfClustering_dS = dataIncompressible.SelfClustering_dS,
       Vs = dataIncompressible.Vs);*/
  //   end substrateFromIncompressible;
  /*
  redeclare function extends substrateFromIdealGas
    "Cast substrate definition from IdealGas (but the state remains in IdealGas) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIdealGas.MolarWeight,
       z = dataIdealGas.z,
       DfG = dataIdealGas.DfG,
       DfH = dataIdealGas.DfH,
       gamma = 1,
       Cp = dataIdealGas.Cp,
       SelfClustering = 0,
       SelfClustering_dH = 0,
       SelfClustering_dS = 0,
       Vs = 1);
       //TODO check specific volume - Vs
  end substrateFromIdealGas;

  redeclare function extends substrateFromIdealGasMSL
    "Cast substrate definition from IdealGasMSL (but the phase remains IdealGasMSL) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIdealGasMSL.data.MM,
       z = dataIdealGasMSL.z,
       DfG = Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(dataIdealGasMSL) - T0*Chemical.Interfaces.IdealGasMSL.molarEntropyPure(dataIdealGasMSL),
       DfH = Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(dataIdealGasMSL),
       gamma = 1,
       Cp = Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp(dataIdealGasMSL),
       SelfClustering = 0,
       SelfClustering_dH = 0,
       SelfClustering_dS = 0,
       Vs = 1);
       //TODO check specific volume - Vs
  end substrateFromIdealGasMSL;

  redeclare function extends substrateFromIdealGasShomate
    "Cast substrate definition from IdealGasShomate (but the phase remains IdealGasShomate) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIdealGasShomate.MolarWeight,
       z = dataIdealGasShomate.z,
       DfG = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(dataIdealGasShomate) - T0*Chemical.Interfaces.IdealGasShomate.molarEntropyPure(dataIdealGasShomate),
       DfH = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(dataIdealGasShomate),
       gamma = 1,
       Cp = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(dataIdealGasShomate),
       SelfClustering = 0,
       SelfClustering_dH = 0,
       SelfClustering_dS = 0,
       Vs = 1);
       //TODO check specific volume - Vs
       end 
       IdealGasShomate;

*/

    redeclare function extends activityCoefficient
      "Return activity coefficient of the substance in the solution"
    algorithm
      activityCoefficient := substanceData.gamma;
    end activityCoefficient;

    redeclare function extends chargeNumberOfIon
      "Return charge number of the substance in the solution"
    algorithm
      chargeNumberOfIon := substanceData.z;
    end chargeNumberOfIon;

    redeclare function extends molarEnthalpyElectroneutral
      "Molar enthalpy of the pure electroneutral substance"
    algorithm
      //Molar enthalpy:
      // - temperature and pressure shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - dp*Vm)
      //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - (dp/dT)*Vm. As a result dH = dT*Cv - dp*Vm for incompressible substances.

      molarEnthalpyElectroneutral := substanceData.DfH + (T - 298.15)*
        substanceData.Cp;
      //   - (p - 100000) * molarVolumePure(substanceData,T,p,v,I);
    end molarEnthalpyElectroneutral;

    redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
    algorithm
      //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
      //+ substanceData.Cv*log(T/298.15);

      //Molar entropy shift:
      // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
      // - pressure shift: with constant molar volume at constant temperature Vm*dP = -T*dS (small amount of work)
      molarEntropyPure := substanceData.Cp*log(T/298.15)
                                + ((substanceData.DfH - substanceData.DfG)/298.15)
                                                         - (molarVolumePure(
          substanceData,
          T,
          p,
          v,
          I)/T)*(p - 100000);

      //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
      //As data: http://www1.lsbu.ac.uk/water/water_phase_diagram.html
      //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
    end molarEntropyPure;

    redeclare function molarVolumePure
      "Molar volume of the pure substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.MolarVolume molarVolumePure "Molar volume";
    algorithm
      molarVolumePure := substanceData.MolarWeight*substanceData.Vs;
      //incompressible
    end molarVolumePure;

    redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
    algorithm
      molarHeatCapacityCp := substanceData.Cp;
    end molarHeatCapacityCp;

    redeclare function extends molarMassOfBaseMolecule
      "Molar mass of the substance"
    algorithm
      molarMass := substanceData.MolarWeight;
    end molarMassOfBaseMolecule;

    redeclare function selfClustering
      "returns true if substance molecules are joining together to clusters"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Boolean selfClustering;
    algorithm
      selfClustering := (substanceData.SelfClustering>0.5);
    end selfClustering;

    redeclare function selfClusteringBondEnthalpy
      "Enthalpy of joining two base molecules of the substance together to cluster"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEnthalpy selfClusteringEnthalpy;
    algorithm
      selfClusteringEnthalpy := substanceData.SelfClustering_dH;
    end selfClusteringBondEnthalpy;

    redeclare function selfClusteringBondEntropy
      "Entropy of joining two base molecules of the substance together to cluster"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      output Modelica.Units.SI.MolarEntropy selfClusteringEntropy;
    algorithm
      selfClusteringEntropy := substanceData.SelfClustering_dS;
    end selfClusteringBondEntropy;

    redeclare replaceable function specificAmountOfParticles
    "Amount of substance particles per its mass"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfSubstance(unit="mol/kg") "Amount of substance particles per its mass";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
      Real amountOfBaseMolecules,amountOfFreeMolecule,amountOfParticles;
      Real x;
    algorithm
      if not selfClustering(substanceData) then
        specificAmountOfSubstance := 1/substanceData.MolarWeight;
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

        amountOfBaseMolecules:=mass/substanceData.MolarWeight;
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        amountOfFreeMolecule := x*nSolution;

        amountOfParticles := amountOfFreeMolecule/(1 - SelfClustering_K*x);

        specificAmountOfSubstance := amountOfParticles/mass;

        //specificAmountOfSubstance := 1/((SelfClustering_K + 1)*substanceData.MolarWeight);
      end if;
    end specificAmountOfParticles;

    redeclare function specificAmountOfFreeBaseMolecule
      "Amount of substance free base molecule per mass of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0 "Ionic strengh (mole fraction based)";
      input Modelica.Units.SI.Mass mass=1 "Mass of substance in solution";
      input Modelica.Units.SI.AmountOfSubstance nSolution=1 "Amount of substances in solution";
      output Real specificAmountOfFreeBaseMolecule(unit="mol/kg")
        "Amount of substance free base molecule per substance mass";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K,amountOfBaseMolecules,x;
    algorithm
      if not selfClustering(substanceData) then
        specificAmountOfFreeBaseMolecule := 1/substanceData.MolarWeight;
      else
        SelfClustering_dG :=selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);

        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));

        amountOfBaseMolecules:=mass/substanceData.MolarWeight;
        x:=((2*SelfClustering_K+nSolution/amountOfBaseMolecules) -
         sqrt((4*SelfClustering_K*nSolution/amountOfBaseMolecules)+
         (nSolution/amountOfBaseMolecules)^2)) / (2*(SelfClustering_K^2));

        specificAmountOfFreeBaseMolecule := (x*nSolution)/mass;

      end if;
      annotation (Inline=true, smoothOrder=2);
    end specificAmountOfFreeBaseMolecule;

    redeclare replaceable function specificEnthalpy
      "Specific molar enthalpy of the substance with electric potential dependence"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificEnthalpy specificEnthalpy
        "Specific enthalpy";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificEnthalpy := (molarEnthalpy(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondEnthalpy(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

      annotation (Inline=true, smoothOrder=2);
    end specificEnthalpy;

    redeclare replaceable function specificVolume
      "Specific volume of the substance"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";

      output Modelica.Units.SI.SpecificVolume specificVolume
        "Specific volume";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificVolume := (molarVolume(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondVolume(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);
    end specificVolume;

    redeclare replaceable function specificHeatCapacityCp
      "Specific heat capacity at constant pressure"
      extends Modelica.Icons.Function;
      input SubstanceData substanceData "Data record of substance";
      input Modelica.Units.SI.Temperature T=298.15 "Temperature";
      input Modelica.Units.SI.Pressure p=100000 "Pressure";
      input Modelica.Units.SI.ElectricPotential v=0
        "Electric potential of the substance";
      input Modelica.Units.SI.MoleFraction I=0
        "Ionic strengh (mole fraction based)";
      output Modelica.Units.SI.SpecificHeatCapacity specificHeatCapacityCp
        "Specific heat capacity at constant pressure";
    protected
      Modelica.Units.SI.MolarEnergy SelfClustering_dG;
      Real SelfClustering_K;
    algorithm
      if selfClustering(substanceData) then
        SelfClustering_dG := selfClusteringBondEnthalpy(substanceData) - T*
          selfClusteringBondEntropy(substanceData);
        SelfClustering_K := exp(-SelfClustering_dG/(Modelica.Constants.R*T));
      end if;

      specificHeatCapacityCp := (molarHeatCapacityCp(
          substanceData,
          T,
          p,
          v,
          I) + (if selfClustering(substanceData) then
        selfClusteringBondHeatCapacityCp(substanceData)*SelfClustering_K/(
        SelfClustering_K + 1) else 0))/molarMassOfBaseMolecule(substanceData);

        //TODO: + selfClusteringBondEnthalpy * der(K/(K + 1))/der(T) .. if (selfClusteringBondHeatCapacityCp!=0)
    end specificHeatCapacityCp;

    redeclare function extends temperature
      "Temperature of substance from its enthalpy"
    protected
      Modelica.Units.SI.SpecificEnthalpy baseSpecificEnthalpy;
    algorithm

      baseSpecificEnthalpy := specificEnthalpy(
          substanceData,
          298.15,
          p,
          v,
          I);

      T := 298.15 + (h - baseSpecificEnthalpy)/specificHeatCapacityCp(
        substanceData);
    end temperature;

    redeclare function extends solution_temperature
      "Temperature of the solution from enthalpies os substances"
      // Modelica.Units.SI.MoleFraction x[size(X, 1)];
    protected
      Modelica.Units.SI.SpecificEnthalpy solution_h_base;
    /*  Modelica.Units.SI.SpecificHeatCapacity solution_Cp=sum(X[i]*
      substanceData[i].Cp/molarMassOfBaseMolecule(substanceData[i]) for
      i in 1:size(X, 1));*/
    algorithm
      solution_h_base := X*specificEnthalpy(
          substanceData,
          298.15,
          p,
          v,
          I);
      T := 298.15 + (h - solution_h_base)/(X*specificHeatCapacityCp(substanceData));
    end solution_temperature;

     redeclare function extends density
      "Return density of the substance in the solution"
     algorithm
      density := 1/substanceData.Vs;
     end density;

    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end Incompressible;

  package IdealGas "Ideal gas with constant heat capacity"
     extends StateOfMatter;

    constant String StateName="IdealGas"; //type checking
    constant Modelica.Units.SI.Temperature T0=298.15;
    constant Modelica.Units.SI.Pressure p0=100000;

     redeclare record extends SubstanceDataParameters "Base substance data"

      parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=1 "Molar weight of the substance";

      parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")= 0
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")= 0
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      parameter Modelica.Units.SI.ActivityCoefficient gamma=1
      "Activity coefficient of the substance";

      parameter Modelica.Units.SI.MolarHeatCapacity Cp=1
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";


      annotation ( preferredView = "info", Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
     end SubstanceDataParameters;

     redeclare record extends SubstanceData "Base substance data"

      Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa") "Molar weight of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

      Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

      Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

      Modelica.Units.SI.MolarHeatCapacity Cp
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";

       encapsulated operator 'constructor'
             import Modelica.Units.SI.MolarMass;
             import SD = Chemical.Interfaces.IdealGas.SubstanceData;
             import IncompressibleSD = Chemical.Interfaces.Incompressible.SubstanceData;
             import IdealGasMSLSD = Chemical.Interfaces.IdealGasMSL.SubstanceData;
             import IdealGasShomateSD = Chemical.Interfaces.IdealGasShomate.SubstanceData;
             import Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral;
             import Chemical.Interfaces.IdealGasMSL.molarEntropyPure;
             import Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp;
             import Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral;
             import Chemical.Interfaces.IdealGasShomate.molarEntropyPure;
             import Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp;

             function fromReal
               input Real MolarWeight=1 "Molar weight of the substance";
               input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
               input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";
               input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";
               input Real gamma=1 "Activity coefficient of the substance";
               input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
               output SD data=SD(
                         MolarWeight=MolarWeight,
                         z=z,
                         DfG=DfG,
                         DfH=DfH,
                         gamma=gamma,
                         Cp=Cp);
             algorithm
               annotation (Inline=true);
             end fromReal;

             function fromIncompressible
               input IncompressibleSD dataIncompressible;
               output SD data=SD(
                         MolarWeight=dataIncompressible.MolarWeight,
                         z=dataIncompressible.z,
                         DfG=dataIncompressible.DfG,
                         DfH=dataIncompressible.DfH,
                         gamma=dataIncompressible.gamma,
                         Cp=dataIncompressible.Cp);
             algorithm
               annotation (Inline=true);
             end fromIncompressible;

             function fromIdealGas
               input SD dataIdealGas;
               output SD data=dataIdealGas;
             algorithm
               annotation (Inline=true);
             end fromIdealGas;

             function fromIdealGasMSL
               input IdealGasMSLSD dataIdealGasMSL;
               output SD data=SD(
                         MolarWeight=dataIdealGasMSL.data.MM,
                         z=dataIdealGasMSL.z,
                         DfG=Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(
                           dataIdealGasMSL,
                           298.15,
                           100000,
                           0,
                           0) - 298.15*Chemical.Interfaces.IdealGasMSL.molarEntropyPure(
                           dataIdealGasMSL,
                           298.15,
                           100000,
                           0,
                           0),
                         DfH=Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(
                           dataIdealGasMSL,
                           298.15,
                           100000,
                           0,
                           0),
                         gamma=1,
                         Cp=Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp(
                           dataIdealGasMSL,
                           298.15,
                           100000,
                           0,
                           0));
             algorithm
               annotation (Inline=true);
             end fromIdealGasMSL;


             function fromIdealGasShomate
               input IdealGasShomateSD dataIdealGasShomate;
               output SD data=SD(
                         MolarWeight=dataIdealGasShomate.MolarWeight,
                         z=dataIdealGasShomate.z,
                         DfG=Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
                           dataIdealGasShomate,
                           298.15,
                           100000,
                           0,
                           0) - 298.15*Chemical.Interfaces.IdealGasShomate.molarEntropyPure(
                           dataIdealGasShomate,
                           298.15,
                           100000,
                           0,
                           0),
                         DfH=Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(
                           dataIdealGasShomate,
                           298.15,
                           100000,
                           0,
                           0),
                         gamma=1,
                         Cp=Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(
                           dataIdealGasShomate,
                           298.15,
                           100000,
                           0,
                           0));
             algorithm
               annotation (Inline=true);
             end fromIdealGasShomate;


       end 'constructor';
      annotation ( preferredView = "info", Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
     end SubstanceData;

     redeclare connector InputSubstanceData = input SubstanceData;
     redeclare connector OutputSubstanceData = output SubstanceData;

     redeclare replaceable model extends BaseProperties "Base properties of incompressible substance"
     end BaseProperties;

     redeclare function extends processData
     "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
     algorithm
      processData :=  SubstanceDataParameters(
        MolarWeight = 0,
        z = 0,
        DfG = -Modelica.Constants.R*T0*log(K),
        DfH = dH,
        gamma = 1,
        Cp = dCp);
     end processData;

    redeclare function extends firstProductDefinition
       "Return formation definition of the first product of chemical process"
    algorithm
       firstProductDefinition :=  SubstanceDataParameters(
         MolarWeight = ((p[2:end]*productsData.MolarWeight) - (s*substratesData.MolarWeight))/p[1],
         z = ((p[2:end]*productsData.z) - (s*substratesData.z))/p[1],
         DfG = ((p[2:end]*productsData.DfG) - (s*substratesData.DfG) - processData.DfG)/p[1],
         DfH = ((p[2:end]*productsData.DfH) - (s*substratesData.DfH) - processData.DfH)/p[1],
         gamma = 1,
         Cp = ((p[2:end]*productsData.Cp) - (s*substratesData.Cp) - processData.Cp)/p[1]);
    end firstProductDefinition;

  /*
  redeclare function extends substrateFromIncompressible
    "Cast substrate definition from Incompressible (but the phase remains Incompressinble) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIncompressible.MolarWeight,
       z = dataIncompressible.z,
       DfG = dataIncompressible.DfG,
       DfH = dataIncompressible.DfH,
       gamma = 1,
       Cp = dataIncompressible.Cp);
  end substrateFromIncompressible;

  redeclare function extends substrateFromIdealGas
    "Cast substrate definition from IdealGas (but the state remains in IdealGas) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data := dataIdealGas;
  end substrateFromIdealGas;

  redeclare function extends substrateFromIdealGasMSL
    "Cast substrate definition from IdealGasMSL (but the phase remains IdealGasMSL) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIdealGasMSL.data.MM,
       z = dataIdealGasMSL.z,
       DfG = Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(dataIdealGasMSL) - T0*Chemical.Interfaces.IdealGasMSL.molarEntropyPure(dataIdealGasMSL),
       DfH = Chemical.Interfaces.IdealGasMSL.molarEnthalpyElectroneutral(dataIdealGasMSL),
       gamma = 1,
       Cp = Chemical.Interfaces.IdealGasMSL.molarHeatCapacityCp(dataIdealGasMSL));
  end substrateFromIdealGasMSL;

  redeclare function extends substrateFromIdealGasShomate
    "Cast substrate definition from IdealGasShomate (but the phase remains IdealGasShomate) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data :=  SubstanceData(
       MolarWeight = dataIdealGasShomate.MolarWeight,
       z = dataIdealGasShomate.z,
       DfG = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(dataIdealGasShomate) - T0*Chemical.Interfaces.IdealGasShomate.molarEntropyPure(dataIdealGasShomate),
       DfH = Chemical.Interfaces.IdealGasShomate.molarEnthalpyElectroneutral(dataIdealGasShomate),
       gamma = 1,
       Cp = Chemical.Interfaces.IdealGasShomate.molarHeatCapacityCp(dataIdealGasShomate));
  end substrateFromIdealGasShomate;

*/


   redeclare function extends activityCoefficient
    "Return activity coefficient of the substance in the solution"
   algorithm
       activityCoefficient := substanceData.gamma;
   end activityCoefficient;

   redeclare function extends chargeNumberOfIon
    "Return charge number of the substance in the solution"
   algorithm
      chargeNumberOfIon := substanceData.z;
   end chargeNumberOfIon;

   redeclare function extends molarEnthalpyElectroneutral
    "Molar enthalpy of the pure substance in electroneutral solution"
   algorithm
       //Molar enthalpy:
       // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
       //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
       molarEnthalpyElectroneutral := substanceData.DfH
         +(T-298.15)*(substanceData.Cp);
   end molarEnthalpyElectroneutral;

   redeclare function extends molarEntropyPure
    "Molar entropy of the pure substance"
   algorithm
     //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
     //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

       //Molar entropy:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
       // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
       molarEntropyPure := (substanceData.Cp)*log(T/298.15) - Modelica.Constants.R*log(p/100000) + ((substanceData.DfH
        - substanceData.DfG)/298.15);

       //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
   end molarEntropyPure;

   redeclare function extends molarVolumePure
    "Molar volume of the pure substance"
   algorithm
       molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
   end molarVolumePure;

   redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
   algorithm
       molarHeatCapacityCp := substanceData.Cp;
   end molarHeatCapacityCp;

   redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
   algorithm
       molarMass := substanceData.MolarWeight;
   end molarMassOfBaseMolecule;

   redeclare function extends temperature "Temperature of substance from its enthalpy"
   algorithm
        T := 298.15 + (h-specificEnthalpy(substanceData,298.15,p,v,I))/specificHeatCapacityCp(substanceData);
   end temperature;

   redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
   algorithm
        T := 298.15 + (h-X*specificEnthalpy(
             substanceData,
             298.15,
             p,
             v,
             I))/(X*specificHeatCapacityCp(substanceData));
   end solution_temperature;

   redeclare function extends density
        "Return density of the substance in the solution"
   algorithm
           density := substanceData.MolarWeight * p / (Modelica.Constants.R * T);
   end density;

    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGas;

  package IdealGasMSL "Ideal gas from Modelica Standard Library 3.2"
    extends StateOfMatter;

    constant String StateName="IdealGasMSL"; //type checking
    constant Modelica.Units.SI.Temperature T0=298.15;
    constant Modelica.Units.SI.Pressure p0=100000;

    redeclare record SubstanceDataParameters

      parameter Modelica.Media.IdealGases.Common.DataRecord data=Modelica.Media.IdealGases.Common.SingleGasesData.N2 "Definition of the substance";

      parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    end SubstanceDataParameters;

    redeclare record SubstanceData

      Modelica.Media.IdealGases.Common.DataRecord data "Definition of the substance";

      Modelica.Units.SI.ChargeNumberOfIon z
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    end SubstanceData;

    redeclare connector InputSubstanceData = input SubstanceData;
    redeclare connector OutputSubstanceData = output SubstanceData;

    redeclare replaceable model extends BaseProperties "Base properties of incompressible substance"
    end BaseProperties;


    redeclare function extends processData
     "Process changes of Gibbs energy, enthalpy, volume and heat capacity (products - reactants)"
    algorithm
      assert(false, "Recalculation of IdealGasMSL from process parameters is not supported!");

    end processData;

    redeclare function extends firstProductDefinition
       "Return formation definition of the first product of chemical process"
    algorithm
       firstProductDefinition :=  SubstanceDataParameters(
        data = Modelica.Media.IdealGases.Common.DataRecord(
           name = "Product defiuned by process",
           MM = ((p[2:end]*productsData.data.MM) - (s*substratesData.data.MM))/p[1],
           Hf = ((p[2:end]*productsData.data.Hf) - (s*substratesData.data.Hf) - processData.data.Hf)/p[1],
           H0 = ((p[2:end]*productsData.data.H0) - (s*substratesData.data.H0) - processData.data.H0)/p[1],
           Tlimit = s*substratesData.data.Tlimit/(s*ones(size(s,1))),
           alow = {((p[2:end]*productsData.data.alow[i]) - (s*substratesData.data.alow[i]) - processData.data.alow[i])/p[1] for i in 1:7},
           blow = {((p[2:end]*productsData.data.blow[i]) - (s*substratesData.data.blow[i]) - processData.data.blow[i])/p[1] for i in 1:2},
           ahigh = {((p[2:end]*productsData.data.ahigh[i]) - (s*substratesData.data.ahigh[i]) - processData.data.ahigh[i])/p[1] for i in 1:7},
           bhigh = {((p[2:end]*productsData.data.bhigh[i]) - (s*substratesData.data.bhigh[i]) - processData.data.bhigh[i])/p[1] for i in 1:2},
           R_s = s*substratesData.data.R_s/(s*ones(size(s,1)))),
        z = ((p[2:end]*productsData.z) - (s*substratesData.z))/p[1]);
    end firstProductDefinition;

    /*
  redeclare function extends substrateFromIncompressible
    "Cast substrate definition from Incompressible (but the phase remains Incompressinble) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     assert(false,"Chemical.Interfaces.IdealMSL.substrateFromIncompressible is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
  end substrateFromIncompressible;

  redeclare function extends substrateFromIdealGas
    "Cast substrate definition from IdealGas (but the state remains in IdealGas) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     assert(false,"Chemical.Interfaces.IdealMSL.substrateFromIdealGas is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
  end substrateFromIdealGas;

  redeclare function extends substrateFromIdealGasMSL
    "Cast substrate definition from IdealGasMSL (but the phase remains IdealGasMSL) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     data := dataIdealGasMSL;
  end substrateFromIdealGasMSL;

  redeclare function extends substrateFromIdealGasShomate
    "Cast substrate definition from IdealGasShomate (but the phase remains IdealGasShomate) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
    assert(false,"Chemical.Interfaces.IdealMSL.substrateFromIdealGasShomate is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
  end substrateFromIdealGasShomate;
*/

    redeclare function extends activityCoefficient
      "Return activity coefficient of the substance in the solution"
    algorithm
      activityCoefficient := 1;
      annotation (Inline=true, smoothOrder=2);
    end activityCoefficient;

    redeclare function extends chargeNumberOfIon
      "Return charge number of the substance in the solution"
    algorithm
      chargeNumberOfIon := substanceData.z;
      annotation (Inline=true, smoothOrder=2);
    end chargeNumberOfIon;

    redeclare function extends molarEnthalpyElectroneutral
      "Molar enthalpy of the pure substance in electroneutral solution"
    algorithm
      //Molar enthalpy:
      // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
      //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
      //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
      molarEnthalpyElectroneutral := substanceData.data.MM*
        Modelica.Media.IdealGases.Common.Functions.h_T(
          substanceData.data,
          T,
          false,
          Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C);
      annotation (Inline=true, smoothOrder=2);
    end molarEnthalpyElectroneutral;

    redeclare function extends molarEntropyPure
      "Molar entropy of the pure substance"
    algorithm
      //molarEntropyPure := ((substanceData.data.DfH - substanceData.data.DfG)/298.15)
      //+ (substanceData.data.Cp+Modelica.Constants.R)*log(T/298.15);

      //Molar entropy:
      // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
      // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)

      molarEntropyPure := substanceData.data.MM*(
        Modelica.Media.IdealGases.Common.Functions.s0_T(substanceData.data, T) -
        substanceData.data.R_s*log(p/100000));

      //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
      //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
      annotation (Inline=true, smoothOrder=2);
    end molarEntropyPure;

    redeclare function extends molarVolumePure
      "Molar volume of the pure substance"
    algorithm
      molarVolumePure := substanceData.data.MM*substanceData.data.R_s*T/p;
      //ideal gas
      annotation (Inline=true, smoothOrder=2);
    end molarVolumePure;

    redeclare function extends molarHeatCapacityCp
      "Molar heat capacity of the substance at constant pressure"
    algorithm
      molarHeatCapacityCp := substanceData.data.MM*
        Modelica.Media.IdealGases.Common.Functions.cp_T(substanceData.data, T);
      annotation (Inline=true, smoothOrder=2);
    end molarHeatCapacityCp;

    redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
    algorithm
      molarMass := substanceData.data.MM;
      annotation (Inline=true, smoothOrder=2);
    end molarMassOfBaseMolecule;

    redeclare function extends temperature "Temperature of substance from its enthalpy"
    protected
       function f_nonlinear "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
         extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
         input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
         input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
       algorithm
         y := Modelica.Media.IdealGases.Common.Functions.h_T(data=data, T=u,
         exclEnthForm=false,refChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C)
                - h;
       end f_nonlinear;

    algorithm
       T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
         function f_nonlinear(data=substanceData.data, h=h), 200, 6000);
    end temperature;

    redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
     // Modelica.Units.SI.MolarMass MM=x*substanceData.data.MM
     //   "molar mass of solution";
     // Modelica.Units.SI.MassFraction x_mass[:]=(x .* substanceData.data.MM) ./
     //     MM "mass fractions";
    protected
        Modelica.Media.IdealGases.Common.DataRecord solutionData=
           Modelica.Media.IdealGases.Common.DataRecord(
               name="solution_temperature",
               MM= 1/sum(X./substanceData.data.MM),
               Hf= X*substanceData.data.Hf,
               H0= X*substanceData.data.H0,
               Tlimit = X*substanceData.data.Tlimit,
               alow = X*substanceData.data.alow,
               blow = X*substanceData.data.blow,
               ahigh = X*substanceData.data.ahigh,
               bhigh = X*substanceData.data.bhigh,
               R_s = X*substanceData.data.R_s);
                     //),
            //sum through moles, not masses

    algorithm
        T :=temperature(
            SubstanceDataParameters(data=solutionData, z=X*(substanceData.z ./ molarMassOfBaseMolecule(substanceData))),
            h,
            p,
            v,
            I);
    end solution_temperature;

    redeclare function extends density
      "Return density of the substance in the solution"
    algorithm
      density := p/(substanceData.data.R_s*T);
      annotation (Inline=true, smoothOrder=2);
    end density;
    annotation (Documentation(revisions="<html>
<p><i>2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGasMSL;

  package IdealGasShomate "Ideal gas based on Shomate equations"
     extends StateOfMatter;

   constant String StateName="IdealGasShomate"; //type checking
   constant Modelica.Units.SI.Temperature T0=298.15;
   constant Modelica.Units.SI.Pressure p0=100000;

   redeclare record extends SubstanceDataParameters
     "Base substance data based on Shomate equations http://old.vscht.cz/fch/cz/pomucky/fchab/Shomate.html"

    parameter Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa")=1 "Molar weight of the substance";

    parameter Modelica.Units.SI.ChargeNumberOfIon z=0
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    parameter Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")=0
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    parameter Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")=0
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    parameter Modelica.Units.SI.ActivityCoefficient gamma=1
      "Activity coefficient of the substance";

    parameter Modelica.Units.SI.MolarHeatCapacity Cp=1
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";
   /*
 parameter Boolean SelfClustering = false "Pure substance is making clusters (weak bonds between molecules)";

 parameter Modelica.Units.SI.ChemicalPotential SelfClustering_dH=0
   "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";                                                                    //-20000
 parameter Modelica.Units.SI.MolarEntropy SelfClustering_dS=0
   "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";
*/
        parameter Real B(unit="J.mol-1")=0 "Shomate parameter B";
        parameter Real C(unit="J.mol-1")=0 "Shomate parameter C";
        parameter Real D(unit="J.K.mol-1")=0 "Shomate parameter D";
        parameter Real E(unit="J.K2.mol-1")=0 "Shomate parameter E";
        parameter Real X=0 "Shomate parameter X";
        parameter Real A_(unit="J.K.mol-1")=0 "Shomate parameter A'";
        parameter Real E_(unit="K")=1e-8 "Shomate parameter E'";


      annotation (preferredView = "info", Documentation(revisions="<html>
<p><i>2016-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
   end SubstanceDataParameters;

   redeclare record extends SubstanceData "Base substance data based on Shomate equations http://old.vscht.cz/fch/cz/pomucky/fchab/Shomate.html"

    connector InputSubstanceData = input SubstanceData;

    Modelica.Units.SI.MolarMass MolarWeight(displayUnit="kDa") "Molar weight of the substance";

    Modelica.Units.SI.ChargeNumberOfIon z
      "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

    Modelica.Units.SI.MolarEnergy DfG(displayUnit="kJ/mol")
      "Gibbs energy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    Modelica.Units.SI.MolarEnergy DfH(displayUnit="kJ/mol")
      "Enthalpy of formation of the substance at SATP conditions (25 degC, 1 bar)";

    Modelica.Units.SI.ActivityCoefficient gamma
      "Activity coefficient of the substance";

    Modelica.Units.SI.MolarHeatCapacity Cp
      "Molar heat capacity of the substance at  SATP conditions (25 degC, 1 bar)";

   /*
 Boolean SelfClustering "Pure substance is making clusters (weak bonds between molecules)";

 Modelica.Units.SI.ChemicalPotential SelfClustering_dH
   "Enthalpy of bond between two molecules of substance at 25degC, 1 bar";                                                                    //-20000
 Modelica.Units.SI.MolarEntropy SelfClustering_dS
   "Entropy of bond between twoo molecules of substance at 25degC, 1 bar";
*/
     Real B(unit="J.mol-1") "Shomate parameter B";
     Real C(unit="J.mol-1") "Shomate parameter C";
     Real D(unit="J.K.mol-1") "Shomate parameter D";
     Real E(unit="J.K2.mol-1") "Shomate parameter E";
     Real X "Shomate parameter X";
     Real A_(unit="J.K.mol-1") "Shomate parameter A'";
     Real E_(unit="K") "Shomate parameter E'";


      annotation (preferredView = "info", Documentation(revisions="<html>
<p><i>2016-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
   end SubstanceData;

   redeclare connector InputSubstanceData = input SubstanceData;
   redeclare connector OutputSubstanceData = output SubstanceData;

   redeclare replaceable model extends BaseProperties "Base properties of incompressible substance"
   end BaseProperties;

   /*
   redeclare function extends substrateFromIncompressible
    "Cast substrate definition from Incompressible (but the phase remains Incompressinble) - use only for firstProductDefinition, where phase shift will be included"
   algorithm 
     assert(false,"Chemical.Interfaces.IdealGasShomate.substrateFromIncompressible is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
   end substrateFromIncompressible;

  redeclare function extends substrateFromIdealGas
    "Cast substrate definition from IdealGas (but the state remains in IdealGas) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
     assert(false,"Chemical.Interfaces.IdealGasShomate.substrateFromIdealGas is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
  end substrateFromIdealGas;

  redeclare function extends substrateFromIdealGasMSL
    "Cast substrate definition from IdealGasMSL (but the phase remains IdealGasMSL) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
    assert(false,"Chemical.Interfaces.IdealGasShomate.substrateFromIdealGasMSL is not implemented yet! If you need this please contect me at marek@matfyz.cz or github.");
  end substrateFromIdealGasMSL;

  redeclare function extends substrateFromIdealGasShomate
    "Cast substrate definition from IdealGasShomate (but the phase remains IdealGasShomate) - use only for firstProductDefinition, where phase shift will be included"
  algorithm 
    data := dataIdealGasShomate;
  end substrateFromIdealGasShomate;
*/


   redeclare function extends activityCoefficient
    "Return activity coefficient of the substance in the solution"
   algorithm
       activityCoefficient := substanceData.gamma;
   end activityCoefficient;

   redeclare function extends chargeNumberOfIon
    "Return charge number of the substance in the solution"
   algorithm
      chargeNumberOfIon := substanceData.z;
   end chargeNumberOfIon;

   redeclare function extends molarEnthalpyElectroneutral
    "Molar enthalpy of the pure substance in electroneutral solution, where der(Hm)=cp*der(T)"
    protected
     parameter Real T0=298.15;
     Real t=T/1000;
     parameter Real A=substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;

     parameter Real H=substanceData.DfH
       - 1000*(substanceData.A_/((-1 + exp((1000*substanceData.E_)/T0))*substanceData.E_)
       - (1000*substanceData.E)/T0 + 0.001*A*T0
       + 5.*10^(-7)*substanceData.B*T0^2 + (1/3)*10^(-9)*substanceData.C*T0^3
       + 2.5*10^(-13)*substanceData.D*T0^4 + (1/1000)^(1.5)/1.5 * T0^1.5 * substanceData.X);

   algorithm
       //Molar enthalpy:
       // - temperature shift: to reach internal energy change by added heat (at constant amount of substance) dU = n*(dH-d(p*Vm)) = n*(dH - R*dT)
       //   where molar heat capacity at constant volume is Cv = dU/(n*dT) = dH/dT - R. As a result dH = dT*(Cv+R) for ideal gas.
       //   And the relation with molar heat capacity at constant pressure as Cp=Cv+R makes dH = dT*Cp.
       molarEnthalpyElectroneutral :=
       H + 1000*(A*t + substanceData.B*t^2/2 + substanceData.C*t^3/3
       + substanceData.D*t^4/4 - substanceData.E/t + substanceData.X*t^1.5/1.5
       + substanceData.A_/substanceData.E_/(exp(substanceData.E_/t) - 1));

   end molarEnthalpyElectroneutral;

   redeclare function extends molarEntropyPure
    "Molar entropy of the pure substance, where der(Sm) = cp*der(T)/T"
    protected
     parameter Modelica.Units.SI.Temperature T0=298.15;
     Real t=T/1000;
     parameter Modelica.Units.SI.MolarEntropy A= substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;

     parameter Modelica.Units.SI.MolarEntropy G= (((substanceData.DfH - substanceData.DfG)/T0)
       + (500000.* substanceData.E)/T0^2
       - (1000*substanceData.A_)/((-1 + exp((1000*substanceData.E_)/T0))*substanceData.E_*T0)
       - 0.001*substanceData.B*T0 - 5*10^(-7) * substanceData.C * T0^2
       - (1/3)*10^(-9)*substanceData.D*T0^3 - sqrt(0.004*T0)* substanceData.X
       + (substanceData.A_*log(1 - exp(-((1000*substanceData.E_)/T0))))/substanceData.E_^2
       - A*log(0.001*T0));

   algorithm
     //molarEntropyPure := ((substanceData.DfH - substanceData.DfG)/298.15)
     //+ (substanceData.Cp+Modelica.Constants.R)*log(T/298.15);

       //Molar entropy:
       // - temperature shift: to reach the definition of heat capacity at constant pressure Cp*dT = T*dS (small amount of added heat energy)
       // - pressure shift: to reach the ideal gas equation at constant temperature Vm*dP = -T*dS (small amount of work)
       molarEntropyPure := G
         + A*log(t) + substanceData.B*t + substanceData.C*t^2/2 + substanceData.D*t^3/3
         - substanceData.E/(2*t^2)
         + 2*substanceData.X*t^0.5 + substanceData.A_/substanceData.E_/t/(exp(substanceData.E_/t) - 1)
         - substanceData.A_/substanceData.E_^2*log(1 - exp(-substanceData.E_/t))
       - Modelica.Constants.R*log(p/100000);

   /*    AA*Log[t] + BB*t + CC*t^2/2 + DD*t^3/3 - EE/(2*t^2) + 2*X*t^0.5 + G +
 AAA/EEE/t/(Exp[EEE/t] - 1) - AAA/EEE^2*Log[1 - Exp[-EEE/t]]

 G + AA*Log[t] + BB*t + CC*t^2/2 + DD*t^3/3 - EE/(2*t^2) + 2*X*t^0.5 +
 AAA/EEE/t/(Exp[EEE/t] - 1) - AAA/EEE^2*Log[1 - Exp[-EEE/t]]
 */

       //For example at triple point of water should be T=273K, p=611.657Pa, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-166 J/mol/K
       //At T=298K, p=1bar, DfH(l)-DfH(g)=44 kJ/mol and S(l)-s(g)=-119 J/mol/K
   end molarEntropyPure;

   redeclare function extends molarVolumePure
    "Molar volume of the pure substance"
   algorithm
       molarVolumePure := Modelica.Constants.R*T/p; //ideal gas
   end molarVolumePure;

   redeclare function extends molarHeatCapacityCp
    "Molar heat capacity of the substance at constant pressure"
    protected
     parameter Real T0=298.15;
     Real t=T/1000;
     parameter Real A= substanceData.Cp
       - ((10^6 * substanceData.A_* exp(1000*substanceData.E_)/T0)) / ((-1 + exp((1000*substanceData.E_)/T0))^2 * T0^2)
       - (10^6 * substanceData.E)/T0^2 - 0.001*substanceData.B*T0 - 10^(-6) * substanceData.C * T0^2
       - 10^(-9) * substanceData.D * T0^3 - sqrt(1/1000)* T0^0.5 * substanceData.X;
   algorithm
       molarHeatCapacityCp := (A + substanceData.B*t + substanceData.C*t^2 +
       substanceData.D*t^3 + substanceData.E/t^2 + substanceData.X*t^0.5 +
       substanceData.A_/t^2*exp(substanceData.E_/t)/(exp(substanceData.E_/t)-1)^2);
   end molarHeatCapacityCp;

   redeclare function extends molarMassOfBaseMolecule "Molar mass of the substance"
   algorithm
       molarMass := substanceData.MolarWeight;
   end molarMassOfBaseMolecule;

   redeclare function extends temperature "Temperature of substance from its enthalpy"
    protected
        function f_nonlinear "Solve molarEnthalpy(data,T) for T with given molar enthalpy"
          extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
        input SubstanceDataParameters data "Ideal gas data";
          input Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
        algorithm
          y := specificEnthalpy(data,u)
                 - h;
        end f_nonlinear;

   algorithm
     T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function f_nonlinear(data=substanceData, h=h), 200, 6000);

   end temperature;

    redeclare function extends solution_temperature
    "Temperature of the solution from enthalpies os substances"
    protected
      Modelica.Units.SI.MolarMass bMM[size(X,1)] = molarMassOfBaseMolecule(substanceData);
      //this is gas, so the self-clustering is not included:
      Modelica.Units.SI.MoleFraction x[size(X,1)]=(X./molarMassOfBaseMolecule(substanceData))/
        sum(X./molarMassOfBaseMolecule(substanceData)) "mole fractions of substances";
      SubstanceDataParameters solutionData=SubstanceDataParameters(
          MolarWeight=sum(x[i]*substanceData[i].MolarWeight for i in 1:size(X, 1)),
          z=sum(x[i]*substanceData[i].z for i in 1:size(X, 1)),
          DfG=sum(x[i]*substanceData[i].DfG for i in 1:size(X, 1)),
          DfH=sum(x[i]*substanceData[i].DfH for i in 1:size(X, 1)),
          gamma=sum(x[i]*substanceData[i].gamma for i in 1:size(X, 1)),
          Cp=sum(x[i]*substanceData[i].Cp for i in 1:size(X, 1)),
          B=sum(x[i]*substanceData[i].B for i in 1:size(X, 1)),
          C=sum(x[i]*substanceData[i].C for i in 1:size(X, 1)),
          D=sum(x[i]*substanceData[i].D for i in 1:size(X, 1)),
          E=sum(x[i]*substanceData[i].E for i in 1:size(X, 1)),
          X=sum(x[i]*substanceData[i].X for i in 1:size(X, 1)),
          A_=sum(x[i]*substanceData[i].A_ for i in 1:size(X, 1)),
          E_=sum(x[i]*substanceData[i].E_ for i in 1:size(X, 1)));          //TODO: gamma,X,E_ are only estimations
    algorithm
      assert(abs(sum(X))<1e-5,"sum(X) must be 1");
      T := temperature(solutionData,h,p,v,I);
    end solution_temperature;

   redeclare function extends density
       "Return density of the substance in the solution"
   algorithm
          density := substanceData.MolarWeight * p / (Modelica.Constants.R * T);
   end density;

    annotation (Documentation(revisions="<html>
<p><i>2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end IdealGasShomate;

  connector SolutionPort
  "Only for connecting the one solution their substances. Please, do not use it in different way."

    //enthalpy
  Modelica.Units.SI.Temperature T "Temperature of the solution";
  flow Modelica.Units.SI.EnthalpyFlowRate dH
    "Internal enthalpy change of the solution";

    //pressure
  Modelica.Units.SI.Pressure p "Pressure of the solution";
  flow Modelica.Units.SI.VolumeFlowRate dV
    "Volume change of the solution";

    //electric port
  Modelica.Units.SI.ElectricPotential v
    "Electric potential in the solution";
  flow Modelica.Units.SI.ElectricCurrent i "Change of electric charge";

    //Extensive properties of the solution:

    // The extensive quantities here have not the real physical flows.
    // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

    //amount of substances
  Modelica.Units.SI.AmountOfSubstance n "Amount of the solution";
  flow Modelica.Units.SI.AmountOfSubstance nj
    "Amount of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //mass of substances
  Modelica.Units.SI.Mass m "Mass of the solution";
  flow Modelica.Units.SI.Mass mj
    "Mass of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //volume of substances
  Modelica.Units.SI.Volume V "Volume of the solution";
  flow Modelica.Units.SI.Volume Vj
    "Volume of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //Gibbs energy of substances
  Modelica.Units.SI.Energy G "Free Gibbs energy of the solution";
  flow Modelica.Units.SI.Energy Gj
    "Gibbs energy of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //electric charge of the substance
  Modelica.Units.SI.ElectricCharge Q "Electric charge of the solution";
  flow Modelica.Units.SI.ElectricCharge Qj
    "Electric charge of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

    //ionic strength of substances
  Modelica.Units.SI.MoleFraction I
    "Mole fraction based ionic strength of the solution";
  flow Modelica.Units.SI.MoleFraction Ij
    "Mole-fraction based ionic strength of the substance (fictive flow to calculate total extensive property in solution as sum from all substances)";

  /*  //suport for structural properties
  replaceable package stateOfMatter = StateOfMatter  constrainedby StateOfMatter
  "Substance model to translate data into substance properties"
     annotation (choicesAllMatching = true);*/

    annotation (
    defaultComponentName="solution",
    Documentation(revisions="<html>
<p><i>2015-2016</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Solution port integrates all substances of the solution:</p>
<p>Such as if there are connected together with electric port, thermal port and with port composed with the amont of substance and molar change of substance.</p>
</html>"),
         Icon(graphics={            Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,0},
            fillColor={127,127,0},
            fillPattern=FillPattern.Solid)}),
    Diagram(graphics={
     Text(extent={{-160,110},{40,50}},   lineColor={127,127,0},    textString = "%name",
          fillColor={127,127,0},
          fillPattern=FillPattern.Solid),
                    Rectangle(
            extent={{-40,40},{40,-40}},
            lineColor={127,127,0},
            fillColor={127,127,0},
            fillPattern=FillPattern.Solid,
            lineThickness=1)}));
  end SolutionPort;

  model Total "Summation of all extensible properties per substance"


    SolutionPort solution
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    parameter Boolean ElectricGround = true
    "Is the solution electric potential equal to zero during simulation (if not useElectricPort)?"
      annotation (HideResult=true);

    Modelica.Blocks.Interfaces.RealInput pressure "pressure"
      annotation (Placement(transformation(extent={{-120,58},{-80,100}})));
    Modelica.Blocks.Interfaces.RealInput temperature "temperature"
      annotation (Placement(transformation(extent={{-120,-22},{-80,20}})));
    Modelica.Blocks.Interfaces.RealOutput volume_der "derivation of volume"
      annotation (Placement(transformation(extent={{100,80},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput enthalpy_der "derivation of enthalpy"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput gibbsEnergy "Gibbs Energy of solution"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    Modelica.Blocks.Interfaces.RealOutput charge
                                             "electric charge"
      annotation (Placement(transformation(extent={{100,-40},{120,-20}})));
    Modelica.Blocks.Interfaces.RealOutput volume
                                             "volume"
      annotation (Placement(transformation(extent={{100,-70},{120,-50}})));
    Modelica.Blocks.Interfaces.RealOutput mass
                                             "mass"
      annotation (Placement(transformation(extent={{100,-100},{120,-80}})));
  equation

    solution.p =pressure;
    solution.T =temperature;

    volume_der + solution.dV = 0;
    enthalpy_der + solution.dH = 0;

     //aliases
    solution.G = gibbsEnergy;
    solution.Q = charge;
    solution.V = volume;
    solution.m = mass;

    //electric current to solution must be represented with some substances (e.g. mass flow of electrones)
    if ElectricGround then
      //Solution connected to ground has zero voltage. However, electric current from the solution can varies.
        solution.v = 0;
    end if;
    if (not ElectricGround) then
      //Electrically isolated solution has not any electric current from/to the solution. However, electric potential can varies.
        solution.i = 0;
    end if;

    //Extensive properties of the solution:

    // The extensive quantities here have not the real physical flows.
    // They hack the Kirchhof's flow equation to be counted as the sum from all connected substances in the solution.

    //amount of substances
    solution.n + solution.nj = 0; //total amount of solution is the sum of amounts of each substance

    //mass of substances
    solution.m + solution.mj = 0; //total mass of solution is the sum masses of each substance

    //Gibs energy
    solution.G + solution.Gj = 0; //total free Gibbs energy of solution is the sum of free Gibbs energies of each substance

    //enthalpy
  //  solution.H + solution.Hj = 0;  //total free enthalpy of solution is the sum of enthalpies of each substance

    //ionic strength (mole fraction based)
    solution.I + solution.Ij = 0; //total ionic strength of solution is the ionic strengths of each substance

    //electric charge
    solution.Q + solution.Qj = 0; //total electric charge of solution is the sum of charges of each substance

    //volume
    solution.V + solution.Vj = 0; //total volume of solution is the sum of volumes of each substance

                                                                                                      annotation (
      Documentation(revisions="<html>
<p>2015-2018 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<h4>amountOfSubstances = &int; MolarFlows</h4>
<h4>mass = &int; massChanges</h4>
<h4>volume = &int; volumeChanges</h4>
<h4>freeEnthalpy = &int; EnthalpyChanges</h4>
<h4>freeEntropy = &int; EntropyChanges</h4>
<h4>freeGibbsEnergy = &int; GibbsEnergyChanges</h4>
<p>Integration of all substances together into one homogenous mixture - the solution.</p>
</html>"));
  end Total;

  partial model PartialSolution
    "Base chemical solution as homogenous mixture of the substances (only pressure and electric potential are not defined)"


    outer Modelica.Fluid.System system "System wide properties";

    parameter Boolean ElectricGround = true
    "Is electric potential equal to zero?"
      annotation (Evaluate=true, choices(checkBox=true), Dialog(group="Environment relationships"));

  Modelica.Units.SI.Temperature temperature "Temperature";

  Modelica.Units.SI.Pressure pressure "Pressure";

  Modelica.Units.SI.Volume volume "Current volume of the solution";

  Modelica.Units.SI.Mass mass(stateSelect=StateSelect.prefer)
    "Current mass of the solution";

    Total total(ElectricGround=ElectricGround)
      annotation (Placement(transformation(extent={{74,-96},{94,-76}})));

  protected
  Modelica.Units.SI.Energy gibbsEnergy
    "Gibbs energy of the solution relative to start of the simulation";

  Modelica.Units.SI.HeatFlowRate heatFromEnvironment
    "External heat flow rate";

  Modelica.Units.SI.ElectricCharge charge
    "Current electric charge of the solution";

  Modelica.Units.SI.HeatFlowRate enthalpy_der "derivative of enthalpy";
  Modelica.Units.SI.VolumeFlowRate volume_der "derivative of volume";

  equation

    heatFromEnvironment = enthalpy_der;

    //total inputs - thermodynamic state
    total.pressure = pressure;
    total.temperature = temperature;

    //total outputs = extensible properties
    enthalpy_der = total.enthalpy_der;
    volume_der = total.volume_der;
    gibbsEnergy = total.gibbsEnergy;
    charge = total.charge;
    volume = total.volume;
    mass = total.mass;

  end PartialSolution;

  partial model PartialSolutionWithHeatPort
    "Chemical solution as homogenous mixture of the substances"

    extends Interfaces.PartialSolution(temperature(start=temperature_start));

  parameter Modelica.Units.SI.Temperature temperature_start=system.T_ambient
    "Initial temperature of the solution"
    annotation (Dialog(group="Initialization"));

    parameter Boolean useThermalPort = false "Is thermal port pressent?"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

    parameter Boolean ConstantTemperature = true
    "Is temperature constant (if not useThermalPort)?"
       annotation (Evaluate=true, choices(checkBox=true), Dialog(enable=not useThermalPort, group="Environment relationships"));

    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(T=temperature,
        Q_flow=heatFromEnvironment) if useThermalPort annotation (Placement(
          transformation(extent={{-70,-90},{-50,-70}}), iconTransformation(
            extent={{-62,-104},{-58,-100}})));

  initial equation
    temperature = temperature_start;
  equation

    //thermal
    if (not useThermalPort) and ConstantTemperature then
      //Ideal thermal exchange between environment and solution to reach constant temperature
      der(temperature) = 0;
    end if;
    if (not useThermalPort) and (not ConstantTemperature) then
      //Thermally isolated without any thermal exchange with environment
      heatFromEnvironment = 0;
    end if;

                                                                                                      annotation (
      Icon(coordinateSystem(
            preserveAspectRatio=false, initialScale=1, extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p>2018 by Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",
        info="<html>
<h4>amountOfSolution = &sum; amountOfSubstances</h4>
<h4>mass = &sum; massOfSubstances</h4>
<h4>volume = &sum; volumeOfSubstances</h4>
<h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
<p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
</html>"));
  end PartialSolutionWithHeatPort;

  partial model ConditionalSubstanceFlow
  "Input of substance molar flow vs. parametric substance molar flow"

    parameter Boolean useSubstanceFlowInput = false
    "=true, if substance flow is provided via input"
    annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
            Dialog(__Dymola_compact=true));

  parameter Modelica.Units.SI.MolarFlowRate SubstanceFlow=0
    "Volumetric flow of Substance if useSubstanceFlowInput=false"
    annotation (HideResult=true, Dialog(enable=not
          useSubstanceFlowInput));

    Modelica.Blocks.Interfaces.RealInput substanceFlow(start=SubstanceFlow, final unit="mol/s")=q
      if useSubstanceFlowInput
         annotation (HideResult=true,
         Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={40,40})));

  Modelica.Units.SI.MolarFlowRate q "Current Substance flow";
  equation
    if not useSubstanceFlowInput then
      q = SubstanceFlow;
    end if;

  end ConditionalSubstanceFlow;

  partial model PartialSolutionSensor

    import Chemical.Utilities.Types.SolutionChoice;

   // parameter SolutionChoice solutionFrom = Chemical.Utilities.Types.SolutionChoice.fromSubstrate "Chemical solution"
    parameter SolutionChoice solutionFrom = Chemical.Utilities.Types.SolutionChoice.fromSubstrate "Chemical solution"
        annotation(HideResult=true, Dialog(group="Conditional inputs"));

    parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
      annotation (HideResult=true, Dialog(enable=(solutionFrom == SolutionChoice.fromParameter)));

    Chemical.Interfaces.SolutionPort solution(
        T=solutionState.T,
        p=solutionState.p,
        v=solutionState.v,
        n=solutionState.n,
        m=solutionState.m,
        V=solutionState.V,
        G=solutionState.G,
        Q=solutionState.Q,
        I=solutionState.I,
        i=0,
        dH=0,
        dV=0,
        nj=0,
        mj=0,
        Vj=0,
        Gj=0,
        Qj=0,
        Ij=0)
          if (solutionFrom == SolutionChoice.fromSolutionPort) "To connect substance with solution, where is pressented"
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));

      Chemical.Interfaces.SolutionState solutionState;

  protected

      Chemical.Interfaces.InputSolutionState inputSubstrateSolution=solutionState if (solutionFrom == Chemical.Utilities.Types.SolutionChoice.fromSubstrate);

  equation

    if (solutionFrom == SolutionChoice.fromParameter) then
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

  partial model ConditionalKinetics
    "Input of kinetics coefficient vs. parametric kinetics coefficient"

    parameter Boolean useForwardRateInput=false
      "= true, if forward rate coefficient is provided via input"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),
        Dialog(group="Chemical kinetics", __Dymola_compact=true));

    parameter Real k_forward(unit="mol/s") = 1 "Forward rate coefficient (mole-fraction based)  if useForwardRateInput=false"
      annotation (HideResult=true, Dialog(group="Chemical kinetics", enable=not useForwardRateInput));

    Modelica.Blocks.Interfaces.RealInput kfInput(start=KC, final unit="mol2.s-1.J-1")=
       kf if useForwardRateInput
       annotation ( HideResult=true, Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40}),
                          iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-60,40})));

    Real kf(final unit="mol/s") "Current forward rate coefficient";

  equation
    if not useForwardRateInput then
      kf = k_forward;
    end if;

  end ConditionalKinetics;

  partial model SISO "Base Model with basic flow eqautions for SISO"

    import Chemical.Utilities.Types.InitializationMethods;

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the flow" annotation (HideResult=true, Dialog(tab="Advanced"));
    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(HideResult=true, Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
      annotation(HideResult=true, Dialog(tab= "Initialization"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0=0 "Initial value for n_flow"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(initN_flow == InitializationMethods.state)));
    parameter Chemical.Utilities.Units.MolarFlowAcceleration n_acceleration_0=0 "Initial value for der(n_flow)"
      annotation (HideResult=true, Dialog(tab="Initialization", enable=(initN_flow == InitializationMethods.derivative)));


    Chemical.Interfaces.Rear rear( state_rearwards(u=u_rear_out, h=h_rear_out))
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));

    Chemical.Interfaces.Fore fore( state_forwards(u=u_fore_out, h=h_fore_out))
      annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));


    Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect)=rear.n_flow;

  protected

    // delta potential computation
    Modelica.Units.SI.ChemicalPotential du_fore;
    Modelica.Units.SI.ChemicalPotential du_rear;

    Modelica.Units.SI.MoleFraction x_rear, x_fore;
  /*  Modelica.Units.SI.Concentration c_in, c_out;
  Modelica.Units.SI.Molality b_in, b_out;
  Modelica.Units.SI.MassFraction X_in, X_out;
*/

    outer Chemical.DropOfCommons dropOfCommons;

    Chemical.Interfaces.InputSubstanceState state_rear_in "Input substance state in forwards direction";
    Chemical.Interfaces.InputSubstanceState state_fore_in "Input substance state in rearwards direction";

    //Chemical.Interfaces.SubstanceState state_out "Input substance state in rearwards direction";

    // input state quantities
    //Modelica.Units.SI.ChemicalPotential u_rear_in=rear.state_forwards.u "Chemical potential of substance entering";
    //Modelica.Units.SI.ChemicalPotential u_fore_in=fore.state_rearwards.u "Chemical potential of substance entering";
    //Modelica.Units.SI.MolarEnthalpy h_rear_in=rear.state_forwards.h "Enthalpy of substance enetering";
    //Modelica.Units.SI.MolarEnthalpy h_fore_in=fore.state_rearwards.h "Enthalpy of substance enetering";

    //outlet state quantities
    Modelica.Units.SI.ChemicalPotential u_rear_out "Chemical potential of substance exiting";
    Modelica.Units.SI.ChemicalPotential u_fore_out "Chemical potential of substance exiting";
    Modelica.Units.SI.MolarEnthalpy h_rear_out "Enthalpy of substance exiting";
    Modelica.Units.SI.MolarEnthalpy h_fore_out "Enthalpy of substance exiting";

    Modelica.Units.SI.ChemicalPotential uPure_substrate "Electro-chemical potential of pure substance entering";
    Modelica.Units.SI.ChemicalPotential uPure_product "Electro-chemical potential of pure substance exiting";

    //Chemical.Utilities.Units.URT duRT_fore, duRT_rear;

  initial equation
    if initN_flow == InitializationMethods.state then
      n_flow = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(n_flow) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(n_flow) = 0;
    end if;

  equation

    connect(rear.state_forwards, state_rear_in);
    connect(fore.state_rearwards, state_fore_in);

    u_fore_out = state_rear_in.u + du_fore;
    u_rear_out = state_fore_in.u + du_rear;

  //  du_fore = rear.state_forwards.u - fore.state_forwards.u;
  //  du_rear = rear.state_rearwards.u - fore.state_rearwards.u;


    h_rear_out = state_fore_in.h;
    h_fore_out = state_rear_in.h;


  //  duRT_fore = (rear.state_forwards.u / (Modelica.Constants.R*rear.solution.T)) - (fore.state_forwards.u / (Modelica.Constants.R*fore.solution.T));
  //  duRT_rear = (rear.state_rearwards.u / (Modelica.Constants.R*rear.solution.T)) - (fore.state_rearwards.u / (Modelica.Constants.R*fore.solution.T));


    x_rear = exp(((rear.state_forwards.u - uPure_substrate)./(Modelica.Constants.R*rear.solution.T)));
    x_fore = exp(((fore.state_rearwards.u - uPure_product)./(Modelica.Constants.R*fore.solution.T)));


    uPure_substrate = Chemical.Interfaces.Properties.electroChemicalPotentialPure(
      rear.definition,
      rear.solution);
    uPure_product = Chemical.Interfaces.Properties.electroChemicalPotentialPure(
      fore.definition,
      fore.solution);



    fore.n_flow + rear.n_flow = 0;
    fore.r = rear.r - der(rear.n_flow) * L;






    annotation (Documentation(info="<html>
<u>Interface class for all components with one fore and one rear port and a massflow without a mass storage between.</u>
<u>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation, as well as the clipping of u_fore to u_min. </u>
<u>If u_fore should be lower the u_min, the remaining potential drop is added on the difference in inertial potential r, basically accelerating or decelerating the massflow. </u>
<u>The component offers different initialization methods for the massflow, as well as several parameters used in the equations above. </u>
<u>The clipping of the massflow can be turned off (this should be done by the modeler as a final modificator while extending to hide this option from the enduser).</u>
</html>"));
  end SISO;

  connector Inlet "Inlet"

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

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    InputSubstanceState state "State of substance";
    InputSolutionState solution "State of solution";
    stateOfMatter.InputSubstanceData definition "Definition of substance";

    annotation (Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{-100,100},{-40,0},{-100,-100},{100,0},{-100,100}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{52,0},{-48,50},{-28,0},{-48,-50},{52,0}},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5,
            lineColor={158,66,200})}),
      Documentation(revisions="<html>
<p><i>2023</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>&Delta;r = -der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. electro-chemical inertia</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg;(T) .. chemical potential of pure substance</p>
<p>DfG(T) .. free Gibbs energy of formation of the substance at current temperature T. </p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
  end Inlet;

  connector Outlet "Outlet providing substance and solution definition"

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

    Modelica.Units.SI.ChemicalPotential r "Inertial Electro-chemical potential";
    flow Modelica.Units.SI.MolarFlowRate n_flow  "Molar change of the substance";

    OutputSubstanceState state "State of substance in solution";
    OutputSolutionState solution "State of solution";
    stateOfMatter.OutputSubstanceData definition "Definition of substance";


    annotation ( Icon(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{100,0},{-100,100},{-40,0},{-100,-100},{100,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true), graphics={
          Polygon(
            points={{50,0},{-50,50},{-30,0},{-50,-50},{50,0}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={158,66,200},
            lineThickness=0.5)}),
      Documentation(revisions="<html>
<p><i>2025</i></p>
<p>Marek Matejak </p>
</html>",   info="<html>

<p>Chemical streams:</p>
<h4>u = û + r</h4>
<h4>&Delta;r = -der(q)*L</h4>
<p>u .. electro-chemical potential</p>
<p>û .. steady-state electro-chemical potential</p>
<p>r .. electro-chemical inertia</p>
<p>q .. molar flow rate</p>
<p>L .. electro-chemical inductance</p>

<p>Definition of electro-chemical potential of the substance:</p>
<h4>u(x,T,v) = u&deg;(T) + R*T*ln(gamma*x) + z*F*v</h4>
<h4>u&deg;(T) = DfG(T) = DfH - T * DfS</h4>
<p>where</p>
<p>x .. mole fraction of the substance in the solution</p>
<p>T .. temperature in Kelvins</p>
<p>v .. eletric potential of the solution</p>
<p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
<p>R .. gas constant</p>
<p>F .. Faraday constant</p>
<p>gamma .. activity coefficient</p>
<p>u&deg;(T) .. chemical potential of pure substance</p>
<p>DfG(T) .. free Gibbs energy of formation of the substance at current temperature T. </p>
<p>DfH .. free enthalpy of formation of the substance</p>
<p>DfS .. free entropy of formation of the substance </p>
<p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
</html>"));
  end Outlet;

  partial model SISOOnedirectional "Base Model with basic flow eqautions for SISO"
    import Chemical;
    import Chemical.Utilities.Types.InitializationMethods;

    replaceable package stateOfMatterIn = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter "Substance model of inlet"
      annotation (Dialog(tab="Advanced"), choices(
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatterIn =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    replaceable package stateOfMatterOut = Interfaces.Incompressible constrainedby
      Interfaces.StateOfMatter
    "Substance model of outlet"
      annotation (Dialog(tab="Advanced"), choices(
        choice(redeclare package stateOfMatterOut =
          Chemical.Interfaces.Incompressible  "Incompressible"),
        choice(redeclare package stateOfMatterOut =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOfMatterOut =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOfMatterOut =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
    parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

    parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance of the molar flow" annotation (Dialog(tab="Advanced"));

    Chemical.Interfaces.Inlet inlet(redeclare package stateOfMatter =
          stateOfMatterIn)                                                           annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
    Chemical.Interfaces.Outlet outlet(redeclare package stateOfMatter =
          stateOfMatterOut)                                                               annotation (Placement(transformation(extent={{80,-20},{120,20}})));

    Modelica.Units.SI.MolarFlowRate n_flow(stateSelect=n_flowStateSelect) = inlet.n_flow
        "Molar flow through component";

    Modelica.Units.SI.MoleFraction x_in, x_out;
    Modelica.Units.SI.Concentration c_in, c_out;
    Modelica.Units.SI.Molality b_in, b_out;
    Modelica.Units.SI.MassFraction X_in, X_out;

  protected
    outer Chemical.DropOfCommons dropOfCommons;

    Modelica.Units.SI.ChemicalPotential uPure_in;
    Modelica.Units.SI.ChemicalPotential u_in=inlet.state.u "Electro-chemical potential of substance entering";
    Modelica.Units.SI.ChemicalPotential u0_in "Electro-chemical potential of pure substance entering";
    Modelica.Units.SI.MolarEnthalpy h_in=inlet.state.h "Enthalpy of substance enetering";

    //outlet state quantities
    Modelica.Units.SI.ChemicalPotential uPure_out;
    Modelica.Units.SI.ChemicalPotential u_out "Electro-chemical potential of substance exiting";
    Modelica.Units.SI.ChemicalPotential u0_out "Electro-chemical potential of pure substance exiting";
    Modelica.Units.SI.MolarEnthalpy h_out "Enthalpy of substance exiting";

    Modelica.Units.SI.ChemicalPotential du "Electro-chemical gradient";
    Chemical.Utilities.Units.URT duRT;

  initial equation
    if initN_flow == InitializationMethods.state then
      n_flow = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(n_flow) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(n_flow) = 0;
    end if;
  equation

    //assert(duRT>=0,"SISO suports only forward flow");

    duRT = u_in/(Modelica.Constants.R*inlet.solution.T) - u_out/(Modelica.Constants.R*outlet.solution.T);
    du = u_in - u_out;

    c_in = x_in * inlet.solution.n/inlet.solution.V;
    c_out = x_out * outlet.solution.n/outlet.solution.V;

    b_in = x_in * inlet.solution.n/inlet.solution.m;
    b_out = x_out * outlet.solution.n/outlet.solution.m;

    x_in = exp((u_in-u0_in)/(Modelica.Constants.R*inlet.solution.T));
    x_out = exp((u_out-u0_out)/(Modelica.Constants.R*outlet.solution.T));

    X_in = b_in / stateOfMatterIn.specificAmountOfParticles(
      inlet.definition,
      inlet.solution.T,
      inlet.solution.p,
      inlet.solution.v,
      inlet.solution.I);
    X_out = b_out / stateOfMatterOut.specificAmountOfParticles(
      outlet.definition,
      outlet.solution.T,
      outlet.solution.p,
      outlet.solution.v,
      outlet.solution.I);

    uPure_in = stateOfMatterIn.electroChemicalPotentialPure(
      inlet.definition,
      inlet.solution.T,
      inlet.solution.p,
      inlet.solution.v,
      inlet.solution.I);
    uPure_out = stateOfMatterOut.electroChemicalPotentialPure(
      outlet.definition,
      outlet.solution.T,
      outlet.solution.p,
      outlet.solution.v,
      outlet.solution.I);

    u0_in = uPure_in;
    u0_out = uPure_out;

    inlet.n_flow + outlet.n_flow = 0;
    outlet.r = inlet.r - der(inlet.n_flow) * L;

    outlet.state.u = u_out;
    outlet.state.h = h_out;

    h_out = h_in;

    annotation (Documentation(revisions="<html>
<p><i>2025</i></p>
<p><i>by </i>Marek Matejak, Ph.D.</p>
</html>",   info="<html>
<p>Interface class for all components with an Inlet and an Outlet and a molarflow without a mass storage between.</p>
<p>This class already implements the equations that are common for such components, namly the conservation of mass, the intertance equation. </p>
</html>"));
  end SISOOnedirectional;

  partial model PartialChangeSolution "Substance between different chemical solutions"
    extends Chemical.Interfaces.SISOOnedirectional;



    Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{92,-52},{112,-32}}),   iconTransformation(extent={{92,-52},{112,-32}})));

  equation

    outlet.definition = inlet.definition;

    outlet.solution.T = solution.T;
    outlet.solution.p = solution.p;
    outlet.solution.v = solution.v;
    outlet.solution.n = solution.n;
    outlet.solution.m = solution.m;
    outlet.solution.V = solution.V;
    outlet.solution.G = solution.G;
    outlet.solution.Q = solution.Q;
    outlet.solution.I = solution.I;

    solution.dH = 0;
    solution.i = 0;
    solution.Qj = 0;
    solution.Ij = 0;
    solution.nj = 0;
    solution.mj = 0;
    solution.Vj = 0;
    solution.Gj = 0;
    solution.dV = 0;


    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
  end PartialChangeSolution;

  partial model PartialChangeState
    extends Chemical.Interfaces.SISOOnedirectional;


    Interfaces.SolutionPort solution annotation (Placement(transformation(extent={{92,-52},{112,-32}}),   iconTransformation(extent={{92,-52},{112,-32}})));

    parameter stateOfMatterOut.SubstanceDataParameters outputSubstanceData
        annotation (choicesAllMatching = true);

  equation

    outlet.definition = outputSubstanceData;

    outlet.solution.T = solution.T;
    outlet.solution.p = solution.p;
    outlet.solution.v = solution.v;
    outlet.solution.n = solution.n;
    outlet.solution.m = solution.m;
    outlet.solution.V = solution.V;
    outlet.solution.G = solution.G;
    outlet.solution.Q = solution.Q;
    outlet.solution.I = solution.I;

    solution.dH = 0;
    solution.i = 0;
    solution.Qj = 0;
    solution.Ij = 0;
    solution.nj = 0;
    solution.mj = 0;
    solution.Vj = 0;
    solution.Gj = 0;
    solution.dV = 0;


    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
  end PartialChangeState;

  partial model PartialGasToLiquid "Gas to liquid change."


    extends Interfaces.PartialChangeState(redeclare package stateOfMatterIn =
          stateIn, final outputSubstanceData=substanceDataOut);

    replaceable package stateIn = Interfaces.IdealGasMSL constrainedby
      Interfaces.StateOfMatter "Substance model of inlet"
      annotation (choices(
        choice(redeclare package stateIn =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateIn =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateIn =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));


     parameter stateOfMatterOut.SubstanceDataParameters substanceDataOut
        annotation (choicesAllMatching = true);

  equation



    annotation (Documentation(revisions="<html>
<p><i>2025</i></p>
<p><i>by </i>Marek Matejak, Ph.D.</p>
</html>",   info="<html>
<p>Change substance state of matter from gaseous to liquid.</p>
</html>"));
  end PartialGasToLiquid;

  partial model PartialLiquidToGas "Gas to liquid change."

    extends Interfaces.PartialChangeState( redeclare package stateOfMatterOut =
          stateOut, final outputSubstanceData=substanceDataOut);



     replaceable package stateOut = Interfaces.IdealGasMSL constrainedby
      Interfaces.StateOfMatter "Substance model of inlet"
      annotation (choices(
        choice(redeclare package stateOut =
          Chemical.Interfaces.IdealGas        "Ideal Gas"),
        choice(redeclare package stateOut =
          Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
        choice(redeclare package stateOut =
          Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    parameter stateOfMatterOut.SubstanceDataParameters substanceDataOut
        annotation (choicesAllMatching = true);

  equation

    annotation (Documentation(revisions="<html>
<p><i>2025</i></p>
<p><i>by </i>Marek Matejak, Ph.D.</p>
</html>",   info="<html>
<p>Change substance state of matter from gaseous to liquid.</p>
</html>"));
  end PartialLiquidToGas;

  partial model PartialProcess "Abstract chemical process"
    import Chemical;
    import Chemical.Utilities.Types.InitializationMethods;

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

    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow =Chemical.Utilities.Types.InitializationMethods.none  "Initialization method for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow"));
    parameter Modelica.Units.SI.MolarFlowRate n_flow_0 = 0 "Initial value for n_flow"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.state)));
    parameter Utilities.Units.MolarFlowAcceleration n_acceleration_0 = 0 "Initial value for der(n_flow)"
      annotation(Dialog(tab= "Initialization", group="Molar flow", enable=(initN_flow == InitializationMethods.derivative)));

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
      annotation(Dialog(tab="Advanced"));

    parameter Integer nS=0 "Number of substrate types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber s[nS]=ones(nS)
      "Stoichiometric reaction coefficient for substrates"
      annotation (HideResult=true);

    parameter Integer nP=0 "Number of product types"
      annotation ( HideResult=true, Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

    parameter Modelica.Units.SI.StoichiometricNumber p[nP]=ones(nP)
      "Stoichiometric reaction coefficients for products"
      annotation (HideResult=true);

    Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Reaction molar flow rate";

    Chemical.Interfaces.Inlet substrates[nS](redeclare package stateOfMatter =
          stateOfMatter)                                                                    annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-100,0}), iconTransformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-100,0})));

    Chemical.Interfaces.Outlet products[nP](redeclare package stateOfMatter =
          stateOfMatter)                                                                     annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={100,0}), iconTransformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={100,0})));

    Modelica.Units.SI.MolarEnthalpy h_mix;

    Real duRT, du, dr, Sx,Px,Kx;

    Modelica.Units.SI.ChemicalPotential uPure_substrates[nS];
    Modelica.Units.SI.ChemicalPotential uPure_products[nP];
  protected
    outer DropOfCommons dropOfCommons;
    //Modelica.Units.SI.ChemicalPotential du;

  initial equation
    if initN_flow == InitializationMethods.state then
      rr = n_flow_0;
    elseif initN_flow == InitializationMethods.derivative then
      der(rr) = n_acceleration_0;
    elseif initN_flow == InitializationMethods.steadyState then
      der(rr) = 0;
    end if;

  equation
    //the main equation
    //assert(duRT>=0,"MIMO suports only forward flow");

    duRT = ((s * (substrates.state.u ./ (Modelica.Constants.R*substrates.solution.T))) - (p * (products.state.u ./ (Modelica.Constants.R*products.solution.T))));
    du = (s * substrates.state.u) - (p * products.state.u);

    for i in 1:nS loop
     uPure_substrates[i] = stateOfMatter.electroChemicalPotentialPure(
      substrates[i].definition,
      substrates[i].solution.T,
      substrates[i].solution.p,
      substrates[i].solution.v,
      substrates[i].solution.I);
    end for;

    for i in 1:nP loop
     uPure_products[i] = stateOfMatter.electroChemicalPotentialPure(
     products[i].definition,
     products[i].solution.T,
     products[i].solution.p,
     products[i].solution.v,
     products[i].solution.I);
    end for;

    Sx = exp(s * ((substrates.state.u - uPure_substrates)./(Modelica.Constants.R*substrates.solution.T)));
    Px = exp((p * ((products.state.u - uPure_products)./(Modelica.Constants.R*products.solution.T))));
    Kx = exp(- ((s * ((uPure_substrates)./(Modelica.Constants.R*substrates.solution.T))) - (p * ((uPure_products)./(Modelica.Constants.R*products.solution.T)))));

    //reaction molar rates
    rr*s = substrates.n_flow;
    rr*p = -products.n_flow;

    products.state.h = h_mix*ones(nP);

    if
      (rr>0) then
      h_mix*(products.n_flow*ones(nP)) + substrates.n_flow*substrates.state.h = 0;
    else
      h_mix = 0;
    end if;

    dr = (s * substrates.r) - (p * products.r);

    if nP>0 then
      (p * products.r) = (s * substrates.r)  -  der(rr)*L;

      for i in 2:nP loop
        //first product is based on inertial potential,
        //other products are provided as source
        der(products[i].state.u).*TC = products[i].r;
      end for;
    end if;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p><i>2013-2025 by </i>Marek Mateják </p>
</html>",   info="<html>
<h4><span style=\"color: #008000\">Notations</span></h4>
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
  end PartialProcess;

  partial model PartialProcessWithSubstanceData "Chemical process with products definitions"
    extends Chemical.Interfaces.PartialProcess;

    parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP]
     annotation (choicesAllMatching = true);

  equation

    products.definition = productsSubstanceData;
    if nS>0 then
      products.solution = fill(substrates[1].solution,nP);
    end if;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}})),
      Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
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
  end PartialProcessWithSubstanceData;

  connector StateInput = input Chemical.Interfaces.SubstanceState "Substance state as connector"
    annotation (
      defaultComponentName="u",
      Icon(graphics={
        Polygon(
          lineColor={162,29,33},
          fillColor={162,29,33},
          fillPattern=FillPattern.Solid,
          points={{-100.0,100.0},{100.0,0.0},{-100.0,-100.0}})},
        coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}},
          preserveAspectRatio=true,
          initialScale=0.2)),
      Diagram(
        coordinateSystem(preserveAspectRatio=true,
          initialScale=0.2,
          extent={{-100.0,-100.0},{100.0,100.0}}),
          graphics={
        Polygon(
          lineColor={162,29,33},
          fillColor={162,29,33},
          fillPattern=FillPattern.Solid,
          points={{0.0,50.0},{100.0,0.0},{0.0,-50.0},{0.0,50.0}}),
        Text(
          textColor={162,29,33},
          extent={{-10.0,60.0},{-10.0,85.0}},
          textString="%name")}),
      Documentation(info="<html>
<p>Connector with one input signal of type Medium.Thermodynamic state. </p>
</html>"));
end Interfaces;
