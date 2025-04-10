within Chemical;
operator record Definition "Definition of a chemical substance or a chemical process"

  Modelica.Media.IdealGases.Common.DataRecord data "Data record of the substance or process";

  Modelica.Units.SI.ChargeNumberOfIon z
  "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";

  encapsulated operator 'constructor'
    import Definition=Chemical.Definition;
    import DataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;

    function fromDataRecord
      input DataRecord data;
      input Real z;
      output Definition result(data=data,z=z);
    algorithm
      annotation (Inline=true);
    end fromDataRecord;

    function fromValues
      input String Name="Substance or process name";
      input Real MM=1 "Molar mass of the substance";
      input Real z=0 "Charge number of the substance (e.g., 0..uncharged, -1..electron, +2..Ca^(2+))";
      input Real DfG=0 "Gibbs energy of formation of the substance at SATP conditions (T0, p0)";
      input Real DfH=0 "Enthalpy of formation of the substance at SATP conditions (T0, p0)";
      input Real gamma=1 "Activity coefficient of the substance";
      input Real Cp=1 "Molar heat capacity of the substance at  SATP conditions (T0, p0)";
      input Real Vs=0.001 "Specific volume of the pure substance (default 1L/kg)";
      input Real T0=298.15 "Base temperature";
      input Real p0=100000 "Base pressure";
      output Definition data=Definition(
                data=DataRecord(
                  name=Name,
                  MM=MM,
                  Hf=DfH/MM,
                  H0=DfH/MM - T0*Cp/MM,
                  Tlimit=1000,
                  alow={0,0,Cp/R,0,0,0,0},
                  blow={DfH/R,(DfH-DfG)/(R*T0)},
                  ahigh={0,0,Cp/R,0,0,0,0},
                  bhigh={DfH/R,(DfH-DfG)/(R*T0)},
                  R_s=R/MM),
                z=z);
    algorithm
      annotation (Inline=true);
    end fromValues;
  end 'constructor';

  encapsulated operator function '+'
    import Definition=Chemical.Definition;
    import DataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;
    input Definition d1;
    input Definition d2;
    output Definition result " = d1 + d2";

  algorithm
    result :=Definition(
        data=DataRecord(
          name="d1+d2",
          MM=d1.data.MM + d2.data.MM,
          Hf=d1.data.Hf + d2.data.Hf,
          H0=d1.data.H0 + d2.data.H0,
          Tlimit=d1.data.Tlimit,
          alow=(((d1.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d1.data.alow)  .+  (((d2.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d2.data.alow),
          blow=(((d1.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d1.data.blow)  .+  (((d2.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d2.data.blow),
          ahigh=(((d1.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d1.data.ahigh)  .+  (((d2.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d2.data.ahigh),
          bhigh=(((d1.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d1.data.bhigh)  .+  (((d2.data.R_s)/(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)) .* d2.data.bhigh),
          R_s=(if (d1.data.MM + d2.data.MM) > 0 then R/(d1.data.MM + d2.data.MM) else R)),
        z=d1.z + d2.z);
        annotation (Inline=true);
  end '+';

  encapsulated operator '-'
    import Definition=Chemical.Definition;
    import DataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;
   function negate
     input Definition d;
     output Definition result " = - d";
   algorithm
     result :=Definition(
        data=DataRecord(
          name="-d",
          MM=-d.data.MM,
          Hf=-d.data.Hf,
          H0=-d.data.H0,
          Tlimit=d.data.Tlimit,
          alow=(((d.data.R_s)/(R/(d.data.MM))) .* d.data.alow),
          blow=(((d.data.R_s)/(R/(d.data.MM))) .* d.data.blow),
          ahigh=(((d.data.R_s)/(R/(d.data.MM))) .* d.data.ahigh),
          bhigh=(((d.data.R_s)/(R/(d.data.MM))) .* d.data.bhigh),
          R_s=R/(-d.data.MM)),
        z= -d.z);
        annotation (Inline=true);
   end negate;

   function substract
    input Definition d1;
    input Definition d2;
    output Definition result " = d1 - d2";
   algorithm
    result :=Definition(
        data=DataRecord(
          name="d1-d2",
          MM=d1.data.MM - d2.data.MM,
          Hf=d1.data.Hf - d2.data.Hf,
          H0=d1.data.H0 - d2.data.H0,
          Tlimit=d1.data.Tlimit,
          alow=(((d1.data.R_s)/(R)) .* d1.data.alow)  .-  (((d2.data.R_s)/(R)) .* d2.data.alow),
          blow=(((d1.data.R_s)/(R)) .* d1.data.blow)  .-  (((d2.data.R_s)/(R)) .* d2.data.blow),
          ahigh=(((d1.data.R_s)/(R)) .* d1.data.ahigh)  .-  (((d2.data.R_s)/(R)) .* d2.data.ahigh),
          bhigh=(((d1.data.R_s)/(R)) .* d1.data.bhigh)  .-  (((d2.data.R_s)/(R)) .* d2.data.bhigh),
          R_s=R),
        z=d1.z - d2.z);
        annotation (Inline=true);
   end substract;
  end '-';

  encapsulated operator '*'
    import Definition=Chemical.Definition;
    import DataRecord=Modelica.Media.IdealGases.Common.DataRecord;
    constant Real R=1.380649e-23*6.02214076e23;

  function scalar
    input Real n=1 "Stoichiometric coefficient";
    input Definition d;
    output Definition result " = n * d";
  algorithm
    result :=Definition(
        data=DataRecord(
          name="n*d",
          MM=n*d.data.MM,
          Hf=n*d.data.Hf,
          H0=n*d.data.H0,
          Tlimit=d.data.Tlimit,
          alow=((n*d.data.R_s)/(R/(n*d.data.MM))) .* d.data.alow,
          blow=((n*d.data.R_s)/(R/(n*d.data.MM))) .* d.data.blow,
          ahigh=((n*d.data.R_s)/(R/(n*d.data.MM))) .* d.data.ahigh,
          bhigh=((n*d.data.R_s)/(R/(n*d.data.MM))) .* d.data.bhigh,
          R_s=R/(n*d.data.MM)),
        z=n*d.z);
        annotation (Inline=true);
  end scalar;

  function vector
    input Real[:] n "Stoichiometric coefficients";
    input Definition[:] d;
    output Definition result " = n * d";
  algorithm
    result :=Definition(
        data=DataRecord(
          name="n*d",
          MM=n*d.data.MM,
          Hf=n*d.data.Hf,
          H0=n*d.data.H0,
          Tlimit=d[1].data.Tlimit,
          alow={sum({(n[i]*(d[i].data.R_s/(if (n*d.data.MM)>0 then R/(n*d.data.MM) else R))) * d[i].data.alow[j] for i in 1:size(n,1)}) for j in 1:7},
          blow={sum({(n[i]*(d[i].data.R_s/(if (n*d.data.MM)>0 then R/(n*d.data.MM) else R))) * d[i].data.blow[j] for i in 1:size(n,1)}) for j in 1:2},
          ahigh={sum({(n[i]*(d[i].data.R_s/(if (n*d.data.MM)>0 then R/(n*d.data.MM) else R))) * d[i].data.ahigh[j] for i in 1:size(n,1)}) for j in 1:7},
          bhigh={sum({(n[i]*(d[i].data.R_s/(if (n*d.data.MM)>0 then R/(n*d.data.MM) else R))) * d[i].data.bhigh[j] for i in 1:size(n,1)}) for j in 1:2},
          R_s=(if (n*d.data.MM)>0 then R/(n*d.data.MM) else R)),
        z=n*d.z);
        annotation (Inline=true);
  end vector;
  end '*';


end Definition;
