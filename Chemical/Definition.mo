within Chemical;
operator record Definition "Definition of a chemical substance or a chemical process"

  Chemical.Interfaces.DataRecord data "Data record of the substance or process";

  encapsulated operator 'constructor'
    import Definition=Chemical.Definition;
    import ModelicaDataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    import DataRecord=Chemical.Interfaces.DataRecord;
    import PhaseType=Chemical.Interfaces.PhaseType;
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
      input Real Vm=if (phase == PhaseType.Gas) then (R*T0)/p0 else 0.001 "Molar volume of the pure substance at SATP conditions (298.15, 1bar) (default 1L/mol)";
      // input Real gamma=1 "Activity coefficient of the substance";
      output Definition result(
                data=DataRecord(
                  MM=MM,
                  Hf=DfH,
                  H0=DfH - T0*Cp,
                  Tlimit=1000,
                  alow={0,0,Cp/R,0,0,0,0},
                  blow={DfH/R,(DfH-DfG)/(R*T0)},
                  ahigh={0,0,Cp/R,0,0,0,0},
                  bhigh={DfH/R,(DfH-DfG)/(R*T0)},
                  z=z,
                  phase=phase,
                  Vm=Vm)
              //    Name=Name,
                  //,gamma=gamma
);
    algorithm
      annotation (Inline=true);
    end fromFormationEnergies;
  end 'constructor';

  encapsulated operator function '+'
    import Definition=Chemical.Definition;
    import DataRecord=Chemical.Interfaces.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;
    input Definition d1;
    input Definition d2;
    output Definition result " = d1 + d2";

  algorithm
    result :=Definition(
        data=DataRecord(
          MM=d1.data.MM + d2.data.MM,
          Hf=d1.data.Hf + d2.data.Hf,
          H0=d1.data.H0 + d2.data.H0,
          Tlimit=d1.data.Tlimit,
          alow=d1.data.alow .+ d2.data.alow,
          blow=d1.data.blow .+ d2.data.blow,
          ahigh=d1.data.ahigh .+ d2.data.ahigh,
          bhigh=d1.data.bhigh .+ d2.data.bhigh,
          z=d1.data.z + d2.data.z,
          phase=d1.data.phase,
          Vm=d1.data.Vm + d2.data.Vm));
         // Name="d1+d2",
        annotation (Inline=true);
  end '+';

  encapsulated operator '-'
    import Definition=Chemical.Definition;
    import DataRecord=Chemical.Interfaces.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;
   function negate
     input Definition d;
     output Definition result " = - d";
   algorithm
     result :=Definition(
        data=DataRecord(
          MM=-d.data.MM,
          Hf=-d.data.Hf,
          H0=-d.data.H0,
          Tlimit=d.data.Tlimit,
          alow=(-1) .* d.data.alow,
          blow=(-1) .* d.data.blow,
          ahigh=(-1) .* d.data.ahigh,
          bhigh=(-1) .* d.data.bhigh,
          z= -d.data.z,
          phase = d.data.phase,
          Vm= -d.data.Vm));
        //  Name="-d",
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
          Tlimit=d1.data.Tlimit,
          alow=d1.data.alow   .-  d2.data.alow,
          blow=d1.data.blow   .-  d2.data.blow,
          ahigh=d1.data.ahigh .-  d2.data.ahigh,
          bhigh=d1.data.bhigh .-  d2.data.bhigh,
          z=d1.data.z - d2.data.z,
          phase=d1.data.phase,
          Vm=d1.data.Vm - d2.data.Vm));
      //    Name="d1-d2",
        annotation (Inline=true);
   end substract;
  end '-';

  encapsulated operator '*'
    import Definition=Chemical.Definition;
    import DataRecord=Chemical.Interfaces.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;

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
          Tlimit=d.data.Tlimit,
          alow=n * d.data.alow,
          blow=n * d.data.blow,
          ahigh=n * d.data.ahigh,
          bhigh=n * d.data.bhigh,
          z=n * d.data.z,
          phase=d.data.phase,
          Vm=n * d.data.Vm));
      //    Name="n*d",
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
          Tlimit=d[1].data.Tlimit,
          alow= {sum({(n[i] * d[i].data.alow[j])  for i in 1:size(n,1)}) for j in 1:7},
          blow= {sum({(n[i] * d[i].data.blow[j])  for i in 1:size(n,1)}) for j in 1:2},
          ahigh={sum({(n[i] * d[i].data.ahigh[j]) for i in 1:size(n,1)}) for j in 1:7},
          bhigh={sum({(n[i] * d[i].data.bhigh[j]) for i in 1:size(n,1)}) for j in 1:2},
          z=n*d.data.z,
          phase=d[1].data.phase,
          Vm=n*d.data.Vm));
     //    Name="n*d",
        annotation (Inline=true);
  end vector;
  end '*';


end Definition;
