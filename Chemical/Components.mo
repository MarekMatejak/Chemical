within Chemical;
package Components "Chemical Components"

  model Substance "Substance in solution"
    extends Icons.Substance;

    Modelica.Units.SI.Concentration c(displayUnit="mmol/l")
      "Molar concentration of particles";

    extends Interfaces.PartialSubstanceInSolutionWithAdditionalPorts;

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
        amountOfBaseMolecules + q*molarEnthalpy - q*actualStream(port_a.h_outflow) + (
        if (calculateClusteringHeat) then stateOfMatter.selfClusteringBondEnthalpy(
        substanceData)*der(amountOfBonds) else 0)
                      "heat transfer from other substances in solution [J/s]";

      solution.Gj =amountOfBaseMolecules*port_a.u + amountOfBonds*SelfClustering_dG
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
                -q*actualStream(port_a.h_outflow)
                "heat transfer from other substances in solution [J/s]";

      solution.Gj = amountOfBaseMolecules*port_a.u "Gibbs energy of the substance [J]";

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
  extends Modelica.Icons.Package;

  model Reaction "Chemical Reaction"
    extends Interfaces.ConditionalKinetics;

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

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient"
      annotation(Dialog(group="Chemical kinetics"));

    Modelica.Units.SI.MolarFlowRate rr(start=0) "Reaction molar flow rate";

    Interfaces.SubstancePorts_b substrates[nS] annotation (Placement(
        transformation(extent={{-10,-40},{10,40}},
        rotation=180,
        origin={-100,0}),                             iconTransformation(
          extent={{-10,-40},{10,40}},
        rotation=180,
        origin={-100,0})));

    Interfaces.SubstancePorts_b products[nP] annotation (Placement(
        transformation(extent={{-10,-40},{10,40}},
        rotation=180,
        origin={100,0}),                            iconTransformation(extent={{-10,-40},
            {10,40}},
        rotation=180,
        origin={100,0})));

    Modelica.Units.SI.MolarEnthalpy h_mix;

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));
  protected
    Modelica.Units.SI.ChemicalPotential du;
  equation
    //the main equation
    du = ((p * products.u) - (s * substrates.u));
    rr = - kC * du * exp(-kE*abs(du));

    //reaction molar rates
    rr*s = substrates.q;
    rr*p = -products.q;

    // Implicit definition of the inStream()operator applied to inside connector i

    substrates.h_outflow = h_mix*ones(nS);
    products.h_outflow = h_mix*ones(nP);

    if
      (rr>0 and not EnthalpyNotUsed) then
      h_mix*(products.q*ones(nP)) + substrates.q*inStream(substrates.h_outflow) = 0;
    elseif
          (rr<0 and not EnthalpyNotUsed) then
      h_mix*(substrates.q*ones(nS)) + products.q*inStream(products.h_outflow) = 0;
    else
      h_mix=0;
    end if;

    // 0 = substrates.q * actualStream(substrates.h_outflow) + products.q * actualStream(products.h_outflow);
  /*  0 = sum(substrates[j].q*(if
                             (substrates[j].q > 0) then h_mix else inStream(substrates[j].h_outflow)) for j in 1:nS)
     +sum(products[k].q * (if
                             (products[k].q > 0)   then h_mix else inStream(products[k].h_outflow)) for k in 1:nP);
*/
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Rectangle(
            extent={{-100,-30},{100,30}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-100,-72},{100,-40}},
            lineColor={128,0,255},
          textString="%name"),
          Polygon(
            points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
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
  end Reaction;

  model ElectronTransfer
  "Electron transfer from the solution to electric circuit"
    extends Icons.ElectronTransfer;
    extends Interfaces.PartialSubstanceInSolution(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible,
      final substanceData = Chemical.Interfaces.Incompressible.SubstanceData(
      MolarWeight=5.4857990946e-7,
      z=-1,
      DfH=0,
      DfG=0,
      Cp=0,
      density=1e20));

    Modelica.Electrical.Analog.Interfaces.PositivePin pin annotation (
        Placement(transformation(extent={{90,50},{110,70}}), iconTransformation(
            extent={{-110,-10},{-90,10}})));

  equation

    //electric
    pin.v = electricPotential;
    pin.i + z*Modelica.Constants.F*port_a.q + solution.i = 0;

    //pure substance
    x = 1;

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

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
          Text(
            extent={{-146,-44},{154,-84}},
            textString="%name",
            lineColor={128,0,255})}),
      Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end ElectronTransfer;

  model Diffusion "Solute diffusion"
    extends Icons.Diffusion;
    extends Interfaces.OnePort;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
  Modelica.Units.SI.ChemicalPotential du;
  equation
    //the main equation
    du = (port_b.u - port_a.u);
    port_b.q = kC * du * exp(-kE*abs(du));

     annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
  end Diffusion;

  model GasSolubility "Henry's law of gas solubility in liquid."

    extends Icons.GasSolubility;

    Interfaces.SubstancePort_b gas_port "Gaseous solution"
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));

    Interfaces.SubstancePort_b liquid_port "Dissolved in liquid solution" annotation (Placement(
          transformation(extent={{-10,-110},{10,-90}}), iconTransformation(extent={{-10,-110},{10,-90}})));

    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J") = 0 "Kinetic turnover coefficient";

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));
  protected
    Modelica.Units.SI.ChemicalPotential du;
  equation
    gas_port.q + liquid_port.q = 0;

    du = (liquid_port.u - gas_port.u);

    liquid_port.q = kC*du*exp(-kE*abs(du));

    gas_port.h_outflow = if (EnthalpyNotUsed) then 0 else inStream(liquid_port.h_outflow);
    liquid_port.h_outflow = if (EnthalpyNotUsed) then 0 else inStream(gas_port.h_outflow);

    annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Gaseuous substance dissolition in liquid (Henry&apos;s law, Raoult&apos;s law, Nernst dissolution in one). </p>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>K<sub>H</sub> =x<sub>L</sub> / x<sub>g</sub>&nbsp;</p></td>
<td><p>Henry&apos;s coefficient, Raoult&apos;s coefficient</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G = &Delta;<sub>f</sub>G<sub>L </sub>- &Delta;<sub>f</sub>G<sub>g </sub>= &Delta;<sub>sol</sub>H - T&middot;&Delta;<sub>sol</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K<sub>H</sub>&middot; (f<sub>L</sub> / f<sub>g</sub>)) </p></td>
<td><p>molar Gibb&apos;s energy of the dissolition</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>H = &Delta;<sub>f</sub>H<sub>L </sub>- &Delta;<sub>f</sub>H<sub>g</sub></p></td>
<td><p>molar enthalpy of the dissolition</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>S = &Delta;<sub>f</sub>S<sub>L</sub> - &Delta;<sub>f</sub>S<sub>g</sub> = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>sol</sub>&omega;) </p></td>
<td><p>molar entropy of the dissolition</p></td>
</tr>
</table>
<h4><span style=\"color:#008000\">Notations</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>x<sub>L</sub></p></td>
<td><p>mole fraction of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>x<sub>g</sub></p></td>
<td><p>mole fraction of the substance in the gas</p></td>
</tr>
<tr>
<td><p>f<sub>L</sub></p></td>
<td><p>activity coefficient of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>f<sub>g</sub></p></td>
<td><p>activity coefficient of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>L</sub></p></td>
<td><p>molar enthalpy of formation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>H<sub>g</sub></p></td>
<td><p>molar enthalpy of formation of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>L</sub></p></td>
<td><p>molar entropy of formation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>f</sub>S<sub>g</sub></p></td>
<td><p>molar entropy of formation of the substance in the gas</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>G</p></td>
<td><p>molar Gibbs energy of dissolvation of the substance in the liquid</p></td>
</tr>
<tr>
<td><p>&Delta;<sub>sol</sub>&omega;</p></td>
<td><p>change of number of microstates of particles by dissolution</p></td>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>
</html>"));
  end GasSolubility;

  model Membrane
  "Passive transport of the substance through semipermeable membrane"
    extends Icons.Membrane;
    extends Interfaces.OnePort;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  protected
  Modelica.Units.SI.ChemicalPotential du;
  equation
    //the main equation
    du = (port_a.u - port_b.u);
    port_a.q = kC * du * exp(-kE*abs(du));

    annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>",
        revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
          Text(
            extent={{-97,-12},{97,12}},
            textString="%name",
            lineColor={128,0,255},
          origin={69,2},
          rotation=90)}));
  end Membrane;

  model SubstancePump "Prescribed sunstance molar flow"
    extends Interfaces.OnePort;
    extends Interfaces.ConditionalSubstanceFlow;

  equation
    port_a.q = q;

   annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics={
          Rectangle(
            extent={{-100,-50},{100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Polygon(
            points={{-80,25},{80,0},{-80,-25},{-80,25}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            rotation=360),
          Text(
            extent={{-150,-20},{150,20}},
            lineColor={128,0,255},
            origin={-10,-76},
            rotation=360,
            textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
  end SubstancePump;

  model Speciation
  "Quaternary macromolecule form defined by all its subunits"
    extends Icons.Speciation;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter Integer NumberOfSubunits=1
    "Number of independent subunits occurring in macromolecule";

    Interfaces.SolutionPort solution                                                              annotation (Placement(transformation(extent={{-70,
              -110},{-50,-90}}),
          iconTransformation(extent={{-70,-110},{-50,-90}})));

    Modelica.Units.SI.AmountOfSubstance nm
      "Amount of the macromolecule (all form in the conformation)";
    Modelica.Units.SI.MoleFraction xm
      "Mole fraction of the macromolecule (all form of in the conformation)";

  public
    Interfaces.SolutionPort subunitSolution "The port to connect all subunits"
      annotation (Placement(transformation(extent={{-70,92},{-50,112}}),
          iconTransformation(extent={{30,50},{50,70}})));
  Interfaces.SubstancePort_a port_a annotation (Placement(transformation(
          extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},
            {110,-90}})));
  Interfaces.SubstancePorts_b subunits[NumberOfSubunits]
    "Subunits of macromolecule" annotation (Placement(transformation(extent={
            {-56,-14},{-36,66}}), iconTransformation(
        extent={{-10,-40},{10,40}},
        rotation=90,
        origin={-30,102})));

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  protected
    Modelica.Units.SI.MolarEnthalpy h_mix;
  equation
    //amount of macromolecule (all forms in conformation)
    nm*NumberOfSubunits + subunitSolution.nj = 0;

    //change of macromolecule = change of its subunits
    subunits.q = -port_a.q * ones(NumberOfSubunits);

    //mole fraction of all forms in conformation
    xm = nm/solution.n;

    //electrochemical potential of the specific form
    port_a.u = Modelica.Constants.R*solution.T*log(xm) +
          sum(subunits.u - Modelica.Constants.R*solution.T*log(xm)
           * ones(NumberOfSubunits));

    port_a.h_outflow = h_mix;
    subunits.h_outflow = (h_mix/NumberOfSubunits)*ones(NumberOfSubunits);

     if
       (port_a.q < 0 and not EnthalpyNotUsed) then
       h_mix = inStream(subunits.h_outflow) * ones(NumberOfSubunits);
     elseif
           (port_a.q > 0 and not EnthalpyNotUsed) then
       h_mix = inStream(port_a.h_outflow);
     else
       h_mix = 0;
     end if;

    //properties from subunits
    subunitSolution.dH + solution.dH = 0;
    subunitSolution.i + solution.i = 0;
    subunitSolution.Qj + solution.Qj = 0;
    subunitSolution.Ij + solution.Ij = 0;

    //properties of macromolecule as a whole
    subunitSolution.nj + solution.nj*NumberOfSubunits = 0; //only amount of substance is necessery to express between sites' solution and real solution
    subunitSolution.mj + solution.mj = 0;
    subunitSolution.Vj + solution.Vj = 0;
    subunitSolution.Gj + solution.Gj = 0;
    subunitSolution.dV + solution.dV = 0;

    //shift global solution status to subunits
    subunitSolution.T = solution.T;
    subunitSolution.v = solution.v;
    subunitSolution.p = solution.p;
    subunitSolution.n = solution.n;
    subunitSolution.m = solution.m;
    subunitSolution.V = solution.V;
    subunitSolution.G = solution.G;
    subunitSolution.Q = solution.Q;
    subunitSolution.I = solution.I;

    annotation (defaultComponentName="macromolecule",
      Documentation(revisions="<html>
<p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic</p>
</html>",   info="<html>
<p><b>Macromolecule speciation in chemical equilibrium</b> </p>
<p>The equilibrium of the conformation reactions of macromolecules can be simplified to the reactions of their selected electro-neutral forms of the selected conformation, because of the law of detailed balance.</p>
<p>The assumptions of this calculation are:</p>
<ol>
<li><i>Initial total concentrations of each subunit must be set to the total macromolecule concentration (for the selected conformation).</i></li>
<li><i>The charge, enthalpy of formation, entropy of formation and molar volume of each selected independent subunit form is zero. </i></li>
<li><i>Subunits are connected to the same solution as the macromolecule. </i></li>
</ol>
<h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>x<sub>m</sub>&nbsp;</p></td>
<td><p>the probability of macromolecule(of the selected conformation)</p></td>
</tr>
<tr>
<td><p>f<sub>i</sub> = (x<sub>i</sub>/x<sub>m</sub>)</p></td>
<td><p>the probalitivy of selected independent subunits forms (of the macromolecule in the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>s </sub>= x<sub>m</sub>&middot; &Pi; f<sub>i</sub> = x<sub>m</sub>&middot; &Pi; (x<sub>i</sub>/x<sub>m</sub>)</p></td>
<td><p>the probability of the selected form of macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>s </sub>= u<sub>s</sub>&deg; + R&middot;T&middot;ln(x<sub>m</sub>) + &sum; (u<sub>i</sub> - R&middot;T&middot;ln(x<sub>m</sub>))</p></td>
<td><p>final equation of the equilibrium of electro-chemical potential</p></td>
</tr>
</table>
<p><br><br><br><b><font style=\"color: #008000; \">Notations</font></b></p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p>n<sub>T</sub></p></td>
<td><p>total amount of substances in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>m</sub></p></td>
<td><p>total amount of the macromolecule (of the selected conformation) in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>s</sub></p></td>
<td><p>amount of the specific form of the macromolecule (of the selected conformation) in the solution</p></td>
</tr>
<tr>
<td><p>n<sub>i</sub></p></td>
<td><p>amount of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
<tr>
<td><p>x<sub>m </sub>= n<sub>m </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of macromolecule (of the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>s </sub>= n<sub>s </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of the selected form of the whole macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>x<sub>i </sub>= n<sub>i </sub>/ n<sub>T</sub></p></td>
<td><p>mole fraction of i-th macromolecule(of the selected conformation)&apos;s subunit form</p></td>
</tr>
<tr>
<td><p>u<sub>s</sub>&deg;</p></td>
<td><p>base chemical potential of the selected form of the macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>s </sub>= u<sub>s</sub>&deg; + R&middot;T&middot;ln(x<sub>s</sub>)</p></td>
<td><p>chemical potential of the selected form of the macromolecule (composed from selected subunits in the selected conformation)</p></td>
</tr>
<tr>
<td><p>u<sub>i</sub>&deg; = 0</p></td>
<td><p>base chemical potential of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
<tr>
<td><p>u<sub>i </sub>= R&middot;T&middot;ln(x<sub>i</sub>)</p></td>
<td><p>chemical potential of the specific form of the i-th macromolecule(of the selected conformation)&apos;s subunit in the solution</p></td>
</tr>
</table>
<p><br><br><br><br>For example: If the macromolecule M has four identical independent subunits and each subunit can occur in two form F1 and F2, then the probability of macromolecule form S composed only from four subunits in form F1 is P(S)=P(M)*P(F1)^4.</p>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}),
          graphics={                                                        Text(
            extent={{-22,-106},{220,-140}},
            lineColor={128,0,255},
            textString="%name")}));
  end Speciation;

  model Stream "Flow of whole solution"
    extends Boundaries.Internal.ConditionalSolutionFlow;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter stateOfMatter.SubstanceData substanceData
    "Definition of the substance"
       annotation (choicesAllMatching = true);

  Interfaces.SubstancePort_b port_b annotation (Placement(transformation(
          extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},
            {-90,10}})));
    Sensors.MoleFractionSensor moleFractionSensor(
       redeclare package stateOfMatter = stateOfMatter,
       substanceData=substanceData)
      annotation (Placement(transformation(extent={{56,-10},{76,10}})));
    Sensors.MoleFractionSensor moleFractionSensor1(
       redeclare package stateOfMatter = stateOfMatter,
       substanceData=substanceData)
      annotation (Placement(transformation(extent={{-56,-10},{-76,10}})));
    SubstancePump substancePump(useSubstanceFlowInput=true,EnthalpyNotUsed=EnthalpyNotUsed)
      annotation (Placement(transformation(extent={{-14,-74},{6,-54}})));
    Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={0,-38})));
    Modelica.Blocks.Logical.GreaterThreshold greaterThreshold annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={0,-4})));
    Modelica.Blocks.Math.Product product
      annotation (Placement(transformation(extent={{-40,-36},{-20,-16}})));
    Modelica.Blocks.Math.Product product1 annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={30,-26})));
  Interfaces.SubstancePort_b port_a
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    Interfaces.SolutionPort solution
      annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));

   parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  equation
    product.u1=q;
    product1.u2=q;
    greaterThreshold.u=q;

    connect(port_b, moleFractionSensor1.port_a) annotation (Line(
        points={{-100,0},{-76,0}},
        color={158,66,200}));
    connect(moleFractionSensor.port_a, port_a) annotation (Line(
        points={{76,0},{100,0}},
        color={158,66,200}));
    connect(moleFractionSensor1.solution, solution) annotation (Line(
        points={{-60,-10},{-60,-100}},
        color={0,128,255}));
    connect(solution, moleFractionSensor.solution) annotation (Line(
        points={{-60,-100},{60,-100},{60,-10}},
        color={0,128,255}));
    connect(substancePump.substanceFlow, switch1.y) annotation (Line(
        points={{0,-60},{0,-49},{-2.22045e-015,-49}},
        color={0,0,127}));
    connect(switch1.u2, greaterThreshold.y) annotation (Line(
        points={{2.22045e-015,-26},{0,-26},{0,-15}},
        color={255,0,255}));
    connect(product1.u1, moleFractionSensor.moleFraction) annotation (Line(
        points={{42,-32},{50,-32},{50,0},{56,0}},
        color={0,0,127}));
    connect(product.u2, moleFractionSensor1.moleFraction) annotation (Line(
        points={{-42,-32},{-50,-32},{-50,0},{-56,0}},
        color={0,0,127}));
    connect(port_b, substancePump.port_a) annotation (Line(
        points={{-100,0},{-86,0},{-86,-64},{-14,-64}},
        color={158,66,200}));
    connect(substancePump.port_b, port_a) annotation (Line(
        points={{6,-64},{84,-64},{84,0},{100,0}},
        color={158,66,200}));
    connect(product.y, switch1.u1) annotation (Line(
        points={{-19,-26},{-8,-26}},
        color={0,0,127}));
    connect(product1.y, switch1.u3) annotation (Line(
        points={{19,-26},{8,-26}},
        color={0,0,127}));
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
  end Stream;
end Components;
