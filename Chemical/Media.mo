within Chemical;
package Media

  package Water_Incompressible "Incompressible water"
    extends Modelica.Media.Interfaces.PartialMedium(
      final mediumName="Water",
      final singleState=true,
      final reducedX=true,
      final fixedX=true,
      reference_T = 310.15,
      reference_p = 101325,
      Temperature(
        min=273,
        max=350,
        start=310.15));

    package stateOfMatter = Chemical.Interfaces.Incompressible
    "Substances model to translate data into substance properties";
    constant Integer nCS=1 "Number of chemical substances";

    constant stateOfMatter.SubstanceData substanceData[nS] = {
      Chemical.Substances.Water_liquid_without_selfClustering()}
       "Definition of the substances";

    redeclare model extends BaseProperties(final standardOrderComponents=true)
      "Base properties of medium"

    equation
      d = stateOfMatter.molarMass(substanceData[1]) / stateOfMatter.molarVolume(substanceData[1],T=T,p=p);
      h = stateOfMatter.molarEnthalpy(substanceData[1],T=T,p=p) / stateOfMatter.molarMass(substanceData[1]);
      u = h - p/d;
      MM = stateOfMatter.molarMass(substanceData[1]);
      R = 8.3144/MM;
      state.p = p;
      state.T = T;
    end BaseProperties;

    redeclare replaceable record ThermodynamicState
      "A selection of variables that uniquely defines the thermodynamic state"
      extends Modelica.Icons.Record;
      AbsolutePressure p "Absolute pressure of medium";
      Temperature T "Temperature of medium";
      annotation (Documentation(info="<html>

</html>"));
    end ThermodynamicState;

    redeclare function extends setState_pTX
      "Return thermodynamic state as function of p, T and composition X or Xi"
    algorithm
      state.p := p;
      state.T := T;
    end setState_pTX;

    redeclare function extends setState_phX
      "Return thermodynamic state as function of p, h and composition X or Xi"
    algorithm
      state.p := p;
      state.T := stateOfMatter.solution_temperature(substanceData,h*stateOfMatter.molarMass(substanceData[1]),{1},p);
    end setState_phX;

    redeclare function extends dynamicViscosity "Return dynamic viscosity"
    algorithm
      eta := (2.414e-5)*10^(247.8/(state.T-140));  //https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
      annotation (Documentation(info="<html>

</html>"));
    end dynamicViscosity;

    redeclare function extends thermalConductivity
      "Return thermal conductivity"
    algorithm
      lambda := 0.6; //google
      annotation (Documentation(info="<html>

</html>"));
    end thermalConductivity;

    redeclare function extends specificEnthalpy "Return specific enthalpy"
    algorithm
      h := stateOfMatter.molarEnthalpy(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
    end specificEnthalpy;

    redeclare function extends specificEntropy "Return specific entropy"
    protected
      Real a "activity of substance";
      Modelica.SIunits.MolarEnergy u "electro-chemical potential of substances in the solution";
    algorithm
      a := stateOfMatter.activityCoefficient(substanceData[1]);

      u := stateOfMatter.chemicalPotentialPure(
          substanceData[1],
          state.T,
          state.p) + Modelica.Constants.R*state.T*log(a);

      s := stateOfMatter.molarEntropy(u,substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
      annotation (Documentation(info="<html>

</html>"));
    end specificEntropy;

    redeclare function extends specificHeatCapacityCp
      "Return specific heat capacity at constant pressure"
    algorithm
      cp := stateOfMatter.molarHeatCapacityCp(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
      annotation (Documentation(info="<html>

</html>"));
    end specificHeatCapacityCp;

    redeclare function extends specificHeatCapacityCv
      "Return specific heat capacity at constant volume"
    algorithm
      cv := stateOfMatter.molarHeatCapacityCv(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
      annotation (Documentation(info="<html>

</html>"));
    end specificHeatCapacityCv;

    redeclare function extends isentropicExponent "Return isentropic exponent"
      extends Modelica.Icons.Function;
    algorithm
      gamma := 23128; //http://twt.mpei.ac.ru/MCS/Worksheets/WSP/WKDiag15.xmcd
      annotation (Documentation(info="<html>

</html>"));
    end isentropicExponent;

    redeclare function extends velocityOfSound "Return velocity of sound"
      extends Modelica.Icons.Function;
    algorithm
      a := 1481; //wikipedia
      annotation (Documentation(info="<html>

</html>"));
    end velocityOfSound;

    redeclare function extends density
    algorithm
      d := stateOfMatter.molarMass(substanceData[1]) / stateOfMatter.molarVolume(substanceData[1]);
    end density;

    redeclare function extends temperature
    algorithm
      T := state.T;
    end temperature;

    redeclare function extends pressure
    algorithm
      p := state.p;
    end pressure;

    function C_outflow
     input Modelica.SIunits.MassFraction x_mass[nCS];
     output Real C_outflow[nC];
    algorithm
      C_outflow := C_default;
      annotation(Inline=true);
    end C_outflow;

    function Xi_outflow
      input Modelica.SIunits.MassFraction x_mass[nCS];
      output Modelica.SIunits.MassFraction Xi[nXi];
    algorithm
      Xi := ones(0);
      annotation(Inline=true);
    end Xi_outflow;

    function x_mass
      input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
      input Real actualStream_C[nC];
      output Modelica.SIunits.MassFraction x_mass[nCS];
    algorithm
      x_mass := {1};
      annotation(Inline=true);
    end x_mass;

    function concentration "Concentration of base substance molecules from Xi and C"
      input ThermodynamicState state;
      input Modelica.SIunits.MassFraction Xi[nXi];
      input Real C[nC];
      output Modelica.SIunits.Concentration concentration[nCS];
    algorithm
      concentration := { 1/stateOfMatter.molarVolume(substanceData[1])};
    end concentration;

    function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
      input Modelica.SIunits.ElectricPotential v=0;
      input Real I=0;
      output SpecificEnthalpy h[nCS];
    algorithm
      h := zeros(nCS);
    end specificEnthalpyOffsets;
    annotation (Documentation(info="<html>
<p>
This package is a <strong>template</strong> for <strong>new medium</strong> models. For a new
medium model just make a copy of this package, remove the
\"partial\" keyword from the package and provide
the information that is requested in the comments of the
Modelica source.
</p>
</html>"));
  end Water_Incompressible;

    package Air_MixtureGasNasa "Air as mixture of O2,CO2,H2O and N2"

      extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
        mediumName="Air",
        data={
          Modelica.Media.IdealGases.Common.SingleGasesData.O2,
          Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
          Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
          Modelica.Media.IdealGases.Common.SingleGasesData.N2},
        fluidConstants={
          Modelica.Media.IdealGases.Common.FluidData.O2,
          Modelica.Media.IdealGases.Common.FluidData.CO2,
          Modelica.Media.IdealGases.Common.FluidData.H2O,
          Modelica.Media.IdealGases.Common.FluidData.N2},
        substanceNames = {"O2", "CO2", "H2O", "N2"},
        excludeEnthalpyOfFormation = false,
        referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
        h_offset = 0,
        reference_X=x_default .* stateOfMatter.molarMass(substanceData) ./ (x_default*stateOfMatter.molarMass(substanceData)),
        T_default=310.15);

      package stateOfMatter = Chemical.Interfaces.IdealGasMSL
      "Substances model to translate data into substance properties";
  protected
      constant Modelica.SIunits.MoleFraction x_default[nCS]=
       {0.21,
        0.0004,
        0.02,
        0.7696} "Initial mole fractions of all substances (Please check: x_default*ones(nCS) = 1)";

  public
      constant Integer nCS=nX "Number of chemical substances";

      constant stateOfMatter.SubstanceData substanceData[nCS] = {
        stateOfMatter.SubstanceData(                data=data[i]) for i in 1:nX}
         "Definition of the substances";

      function C_outflow
       input Modelica.SIunits.MassFraction x_mass[nCS];
       output Real C_outflow[nC];
      algorithm
        C_outflow := C_default;
        annotation(Inline=true);
      end C_outflow;

      function Xi_outflow
        input Modelica.SIunits.MassFraction x_mass[nCS];
        output Modelica.SIunits.MassFraction Xi[nXi];
      algorithm
        Xi := x_mass[1:nXi];
        annotation(Inline=true);
      end Xi_outflow;

      function x_mass
        input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
        input Real actualStream_C[nC];
        output Modelica.SIunits.MassFraction x_mass[nCS];
      algorithm
        x_mass := actualStream_Xi;
        annotation(Inline=true);
      end x_mass;

      function concentration "Concentration of substances from Xi and C"
        input ThermodynamicState state;
        input Modelica.SIunits.MassFraction Xi[nXi];
        input Real C[nC];
        output Modelica.SIunits.Concentration concentration[nCS];
      algorithm
        concentration := Xi*density(state)./stateOfMatter.molarMass(substanceData);
      end concentration;

      function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
        input Modelica.SIunits.ElectricPotential v=0;
        input Real I=0;
        output SpecificEnthalpy h[nCS];
      algorithm
         h := zeros(nCS); //all substances are electroneutral, so chemical and medium enthalpies are equal
      end specificEnthalpyOffsets;
    end Air_MixtureGasNasa;

    package SimpleAir_C

      extends Modelica.Media.Air.SimpleAir(
       extraPropertiesNames={"O2","CO2","H2O","N2"},
       T_default=310.15,
       X_default=ones(nX),
       C_default= x_default ./ MM_default,
       T0 = 298.15);
       //C is amount of substance per total mass of solution in [mol/kg] (not a concentration, neither the molality)
       //formation enthalpies of substances are not included in enthalpy

      package stateOfMatter = Chemical.Interfaces.IdealGas
      "Substances model to translate data into substance properties";

  protected
      constant Modelica.SIunits.MoleFraction x_default[nC]=
       {0.21,
        0.0004,
        0.02,
        0.7696};

      constant Modelica.SIunits.MolarMass MM_default = x_default*Chemical.Interfaces.IdealGas.molarMass(substanceData);

  public
      constant Integer nCS=nC "Number of chemical substances";

      constant stateOfMatter.SubstanceData substanceData[nC] = {
        Chemical.Substances.Oxygen_gas(),
        Chemical.Substances.CarbonDioxide_gas(),
        Chemical.Substances.Water_gas(),
        Chemical.Substances.Nitrogen_gas()}
         "Definition of the substances";

      function C_outflow "Outflow values for extra properties of fluid connector"
        input Modelica.SIunits.MassFraction x_mass[nCS];
        output Real C_outflow[nC];
      algorithm
          C_outflow := x_mass[1:nC] ./ stateOfMatter.molarMass(substanceData);
          annotation(Inline=true);
      end C_outflow;

      function Xi_outflow "Outflow values for mass fracion of fluid connector"
        input Modelica.SIunits.MassFraction x_mass[nCS];
        output Modelica.SIunits.MassFraction Xi[nXi];
      algorithm
        Xi := X_default[1:nXi];
        annotation(Inline=true);
      end Xi_outflow;

      function x_mass "Mass fractions from actual streams of fluid connector"
        input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
        input Real actualStream_C[nC];
        output Modelica.SIunits.MassFraction x_mass[nCS];
      algorithm
        x_mass := actualStream_C .* stateOfMatter.molarMass(substanceData);
        annotation(Inline=true);
      end x_mass;

      function concentration "Concentration of substances from Xi and C"
        input ThermodynamicState state;
        input Modelica.SIunits.MassFraction Xi[nXi];
        input Real C[nC];
        output Modelica.SIunits.Concentration concentration[nCS];
      algorithm
        concentration := C*density(state);
        annotation(Inline=true);
      end concentration;

      function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
        input Modelica.SIunits.ElectricPotential v=0;
        input Real I=0;
        output SpecificEnthalpy h[nCS];
    protected
        parameter SpecificEnthalpy H[nCS] = stateOfMatter.molarEnthalpy(
            substanceData,
            T0,
            100000,
            0,
            0) ./ stateOfMatter.molarMass(substanceData); //enthalpy at T0
      algorithm
        h := H;
      end specificEnthalpyOffsets;

    end SimpleAir_C;

    package SimpleO2Gas_C

      extends Modelica.Media.Air.SimpleAir(
       extraPropertiesNames={"O2"},
       T_default=310.15, X_default=ones(nX), C_default={1/stateOfMatter.molarMass(substanceData[nCS])},
       T0 = 298.15);

      package stateOfMatter = Chemical.Interfaces.IdealGas
      "Substances model to translate data into substance properties";

      constant Integer nCS=nC "Number of chemical substances";

      constant stateOfMatter.SubstanceData substanceData[nC] = {
        Chemical.Substances.Oxygen_gas()}
         "Definition of the substances";

      function C_outflow "Outflow values for extra properties of fluid connector"
        input Modelica.SIunits.MassFraction x_mass[nCS];
        output Real C_outflow[nC];
      algorithm
          C_outflow := x_mass[1:nC] ./ stateOfMatter.molarMass(substanceData);
          annotation(Inline=true);
      end C_outflow;

      function Xi_outflow "Outflow values for mass fracion of fluid connector"
        input Modelica.SIunits.MassFraction x_mass[nCS];
        output Modelica.SIunits.MassFraction Xi[nXi];
      algorithm
        Xi := X_default[1:nXi];
        annotation(Inline=true);
      end Xi_outflow;

      function x_mass "Mass fractions from actual streams of fluid connector"
        input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
        input Real actualStream_C[nC];
        output Modelica.SIunits.MassFraction x_mass[nCS];
      algorithm
        x_mass := actualStream_C .* stateOfMatter.molarMass(substanceData);
        annotation(Inline=true);
      end x_mass;

      function concentration "Concentration of substances from Xi and C"
        input ThermodynamicState state;
        input Modelica.SIunits.MassFraction Xi[nXi];
        input Real C[nC];
        output Modelica.SIunits.Concentration concentration[nCS];
      algorithm
        concentration := C*density(state);
        annotation(Inline=true);
      end concentration;

      function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
        input Modelica.SIunits.ElectricPotential v=0;
        input Real I=0;
        output SpecificEnthalpy h[nCS];
    protected
        parameter SpecificEnthalpy H[nCS] = stateOfMatter.molarEnthalpy(
            substanceData,
            T0,
            100000,
            0,
            0) ./ stateOfMatter.molarMass(substanceData);
      algorithm
        h := H; //enthalpy at T0
      end specificEnthalpyOffsets;
    end SimpleO2Gas_C;

    package EthanolInWater_C

        extends Modelica.Media.Water.StandardWater(
         extraPropertiesNames={"H2O","C2H5OH"},
         singleState=true, T_default=310.15, X_default=ones(nX), C_default=x_default./MM_default);

        package stateOfMatter = Chemical.Interfaces.Incompressible
        "Substances model to translate data into substance properties";

        constant Integer nCS=nC "Number of chemical substances";

  protected
        constant Modelica.SIunits.MoleFraction x_default[nCS]=
        {0.5, 0.5} "initial mole fraction of substances";

        constant Modelica.SIunits.MolarMass MM_default = x_default*stateOfMatter.molarMass(substanceData);

        constant Modelica.SIunits.MolarVolume MV_default = x_default*stateOfMatter.molarVolume(substanceData);

        constant Modelica.SIunits.Density density_default = MM_default / MV_default "initial density";

  public
        constant stateOfMatter.SubstanceData substanceData[nCS] = {
          Substances.Water_liquid(),
          Substances.Ethanol_liquid()}
           "Definition of the substances";

        function C_outflow "Outflow values for extra properties of fluid connector"
          input Modelica.SIunits.MassFraction x_mass[nCS];
          output Real C_outflow[nC];
        algorithm
            C_outflow := x_mass[1:nC] ./ stateOfMatter.molarMass(substanceData);
            annotation(Inline=true);
        end C_outflow;

        function Xi_outflow "Outflow values for mass fracion of fluid connector"
          input Modelica.SIunits.MassFraction x_mass[nCS];
          output Modelica.SIunits.MassFraction Xi[nXi];
        algorithm
          Xi := X_default[1:nXi];
          annotation(Inline=true);
        end Xi_outflow;

        function x_mass "Mass fractions from actual streams of fluid connector"
          input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
          input Real actualStream_C[nC];
          output Modelica.SIunits.MassFraction x_mass[nCS];
        algorithm
          x_mass := actualStream_C .* stateOfMatter.molarMass(substanceData);
          annotation(Inline=true);
        end x_mass;

        function concentration "Concentration of substances from Xi and C"
          input ThermodynamicState state;
          input Modelica.SIunits.MassFraction Xi[nXi];
          input Real C[nC];
          output Modelica.SIunits.Concentration concentration[nCS];
        algorithm
          concentration := C*density(state);
          annotation(Inline=true);
        end concentration;


     function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
        input Modelica.SIunits.ElectricPotential v=0;
        input Real I=0;
        output SpecificEnthalpy h[nCS];
    protected
       parameter ThermodynamicState state = setState_pT(100000,298.15);
       parameter SpecificEnthalpy H[nCS] = stateOfMatter.molarEnthalpy(
            substanceData,
            298.15,
            100000,
            0,
            0) ./ stateOfMatter.molarMass(substanceData) - {specificEnthalpy(state),specificEnthalpy(state)};
     algorithm
        h := H;
     end specificEnthalpyOffsets;
    end EthanolInWater_C;

    package StandardWater_C

        extends Modelica.Media.Water.StandardWater(
         singleState=true, T_default=310.15, X_default=ones(nX));


        package stateOfMatter = Chemical.Interfaces.Incompressible
        "Substances model to translate data into substance properties";
        constant Integer nCS=1 "Number of chemical substances";

  protected
        constant Modelica.SIunits.MoleFraction x_default[nCS]=
        {1} "initial mole fraction of substances";

        constant Modelica.SIunits.MolarMass MM_default = x_default*stateOfMatter.molarMass(substanceData);

        constant Modelica.SIunits.MolarVolume MV_default = x_default*stateOfMatter.molarVolume(substanceData);

        constant Modelica.SIunits.Density density_default = MM_default / MV_default "initial density";

  public
        constant stateOfMatter.SubstanceData substanceData[nCS] = {
          Substances.Water_liquid()}
           "Definition of the substances";

        function C_outflow "Outflow values for extra properties of fluid connector"
          input Modelica.SIunits.MassFraction x_mass[nCS];
          output Real C_outflow[nC];
        algorithm
            C_outflow := C_default;
            annotation(Inline=true);
        end C_outflow;

        function Xi_outflow "Outflow values for mass fracion of fluid connector"
          input Modelica.SIunits.MassFraction x_mass[nCS];
          output Modelica.SIunits.MassFraction Xi[nXi];
        algorithm
          Xi := X_default[1:nXi];
          annotation(Inline=true);
        end Xi_outflow;

        function x_mass "Mass fractions from actual streams of fluid connector"
          input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
          input Real actualStream_C[nC];
          output Modelica.SIunits.MassFraction x_mass[nCS];
        algorithm
          x_mass := {1};
          annotation(Inline=true);
        end x_mass;

        function concentration "Concentration of substances from Xi and C"
          input ThermodynamicState state;
          input Modelica.SIunits.MassFraction Xi[nXi];
          input Real C[nC];
          output Modelica.SIunits.Concentration concentration[nCS];
        algorithm
          concentration := {density(state)}./stateOfMatter.molarMass(substanceData);
          annotation(Inline=true);
        end concentration;

     function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
        input Modelica.SIunits.ElectricPotential v=0;
        input Real I=0;
        output SpecificEnthalpy h[nCS];
    protected
       parameter ThermodynamicState state = setState_pT(100000,298.15);
       parameter SpecificEnthalpy H[nCS] = stateOfMatter.molarEnthalpy(
            substanceData,
            298.15,
            100000,
            0,
            0) ./ stateOfMatter.molarMass(substanceData) - {specificEnthalpy(state),specificEnthalpy(state)};
     algorithm
        h := H;
     end specificEnthalpyOffsets;
    end StandardWater_C;
end Media;
