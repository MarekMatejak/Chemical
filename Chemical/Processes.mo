within Chemical;
package Processes

  model Reaction "Chemical Reaction"
    extends Interfaces.ConditionalKinetics;
    import Chemical.Utilities.Types.InitializationMethods;

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

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient"
      annotation(Dialog(group="Chemical kinetics"));

    Modelica.Units.SI.MolarFlowRate rr(stateSelect=n_flowStateSelect) "Reaction molar flow rate";

    Interfaces.Inlet substrates[nS] annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-100,0}), iconTransformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-100,0})));

    Interfaces.Outlet products[nP] annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={100,0}), iconTransformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={100,0})));

    Modelica.Units.SI.MolarEnthalpy h_mix;

  protected
    outer DropOfCommons dropOfCommons;
    Modelica.Units.SI.ChemicalPotential du;

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
    du = ((p * products.u) - (s * substrates.u));
    rr = - kC * du * exp(-kE*abs(du));

    //reaction molar rates
    rr*s = substrates.n_flow;
    rr*p = -products.n_flow;

    products.h = h_mix*ones(nP);

    if
      (rr>0) then
      h_mix*(products.n_flow*ones(nP)) + substrates.n_flow*substrates.h = 0;
    else
      h_mix = 0;
    end if;

    if nP>0 then
      (p * products.r) = (s * substrates.r)  -  der(rr)*L;
      for i in 2:nP loop
        //first product is based on inertial potential,
        //other products are provided as source
        der(products[i].u).*TC = products[i].r;
      end for;
    end if;

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

  model Process "Electro-chemical process"
    extends Interfaces.SISOFlow;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  equation
    //the main equation

    n_flow = - kC * du * exp(-kE*abs(du));

     annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"),
         Icon(graphics={Rectangle(extent={{-100,40},{100,-40}}, lineColor={28,108,200})}));
  end Process;

  model Diffusion "Solute diffusion"
    extends Icons.Diffusion;
    extends Interfaces.SISOFlow;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  equation
    //the main equation

    n_flow = - kC * du * exp(-kE*abs(du));


     annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
  end Diffusion;

  model GasSolubility "Henry's law of gas solubility in liquid."

    extends Icons.GasSolubility;

    extends Interfaces.SISOFlowVertical;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  equation
    //the main equation

    n_flow = - kC * du * exp(-kE*abs(du));



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

  model Membrane "Passive transport of the substance through semipermeable membrane"
    extends Icons.Membrane;
     extends Interfaces.SISOFlow;
    extends Interfaces.ConditionalKinetics;

    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";

  equation
    //the main equation

    n_flow = - kC * du * exp(-kE*abs(du));


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
    extends Interfaces.SISOFlow;
    extends Interfaces.ConditionalSubstanceFlow;

  equation
    n_flow = q;

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

  model Stream "Flow of whole solution"
    extends Boundaries.Internal.ConditionalSolutionFlow;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter stateOfMatter.SubstanceData substanceData
    "Definition of the substance"
       annotation (choicesAllMatching = true);

  Interfaces.Inlet           inlet  annotation (Placement(transformation(
          extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},
            {-90,10}})));
    Sensors.MoleFractionSensor moleFractionSensor1(
       redeclare package stateOfMatter = stateOfMatter,
       substanceData=substanceData)
      annotation (Placement(transformation(extent={{-56,-10},{-76,10}})));
    SubstancePump substancePump(useSubstanceFlowInput=true) annotation (Placement(transformation(extent={{-20,-72},{0,-52}})));
    Modelica.Blocks.Math.Product product
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-6,-24})));
  Interfaces.Outlet          port_a
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


    connect(inlet, moleFractionSensor1.port_a) annotation (Line(points={{-100,0},{-76,0}}, color={158,66,200}));
    connect(moleFractionSensor1.solution, solution) annotation (Line(
        points={{-60,-10},{-60,-100}},
        color={0,128,255}));
    connect(product.u2, moleFractionSensor1.moleFraction) annotation (Line(
        points={{-12,-12},{-12,0},{-56,0}},
        color={0,0,127}));
    connect(inlet, substancePump.inlet) annotation (Line(
        points={{-100,0},{-82,0},{-82,-62},{-20,-62}},
        color={158,66,200},
        thickness=0.5));
    connect(substancePump.outlet, port_a) annotation (Line(
        points={{0,-62},{86,-62},{86,0},{100,0}},
        color={158,66,200},
        thickness=0.5));
    connect(product.y, substancePump.substanceFlow) annotation (Line(points={{-6,-35},{-6,-58}}, color={0,0,127}));
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

  model SpeciationIn "Quaternary macromolecule form defined by all its subunits"
    extends Icons.Speciation;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter Integer NumberOfSubunits=1
    "Number of independent subunits occurring in macromolecule";

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
      annotation(Dialog(tab="Advanced"));


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
  Interfaces.Inlet inlet annotation (Placement(transformation(
          extent={{110,-110},{90,-90}}), iconTransformation(extent={{110,-110},{90,-90}})));
  Interfaces.Outlet subunits[NumberOfSubunits]
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
    outer DropOfCommons dropOfCommons;
    Modelica.Units.SI.MolarEnthalpy h_out;
    Modelica.Units.SI.ChemicalPotential u_out;
  equation

    inlet.r + inlet.u = u_out;

    if NumberOfSubunits>0 then
      (ones(NumberOfSubunits) * subunits.r) = (inlet.r)  -  der(inlet.n_flow)*L;
      for i in 2:NumberOfSubunits loop
        //first subunit is based on inertial potential,
        //other subunits are provided as source
        der(subunits[i].u).*TC = subunits[i].r;
      end for;
    end if;


    //amount of macromolecule (all forms in conformation)
    nm*NumberOfSubunits + subunitSolution.nj = 0;

    //change of macromolecule = change of its subunits
    subunits.n_flow = -inlet.n_flow * ones(NumberOfSubunits);

    //mole fraction of all forms in conformation
    xm = nm/solution.n;

    //electrochemical potential of the specific form
    u_out = Modelica.Constants.R*solution.T*log(xm) +
          sum(subunits.u - Modelica.Constants.R*solution.T*log(xm)
           * ones(NumberOfSubunits));

    h_out = inlet.h;
    subunits.h = (inlet.h/NumberOfSubunits)*ones(NumberOfSubunits);


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
  end SpeciationIn;

  model SpeciationOut "Quaternary macromolecule form defined by all its subunits"
    extends Icons.Speciation;

    replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby Interfaces.StateOfMatter
    "Substance model to translate data into substance properties"
       annotation (choicesAllMatching = true);

    parameter Integer NumberOfSubunits=1
    "Number of independent subunits occurring in macromolecule";

    parameter Modelica.Units.SI.Time TC=0.1 "Time constant for electro-chemical potential adaption" annotation (Dialog(tab="Advanced"));
    parameter Utilities.Units.Inertance L = dropOfCommons.L "Inertance of the flow"
      annotation(Dialog(tab="Advanced"));

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
  Interfaces.Outlet outlet annotation (Placement(transformation(
          extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},{110,-90}})));
  Interfaces.Inlet subunits[NumberOfSubunits]
    "Subunits of macromolecule" annotation (Placement(transformation(extent={
            {-56,-14},{-36,66}}), iconTransformation(
        extent={{10,-40},{-10,40}},
        rotation=90,
        origin={-30,102})));

    parameter Boolean EnthalpyNotUsed=false annotation (
      Evaluate=true,
      HideResult=true,
      choices(checkBox=true),
      Dialog(tab="Advanced", group="Performance"));

  protected
    outer DropOfCommons dropOfCommons;
    Modelica.Units.SI.MolarEnthalpy h_out;
    Modelica.Units.SI.ChemicalPotential u_out;
  equation

    outlet.u = u_out;

    if NumberOfSubunits>0 then
      (ones(NumberOfSubunits) * subunits.r) = (outlet.r)  -  der(outlet.n_flow)*L;
    end if;

    //amount of macromolecule (all forms in conformation)
    nm*NumberOfSubunits + subunitSolution.nj = 0;

    //change of macromolecule = change of its subunits
    subunits.n_flow = -outlet.n_flow * ones(NumberOfSubunits);

    //mole fraction of all forms in conformation
    xm = nm/solution.n;

    //electrochemical potential of the specific form
    u_out = Modelica.Constants.R*solution.T*log(xm) +
          sum(subunits.u - Modelica.Constants.R*solution.T*log(xm)
           * ones(NumberOfSubunits));

    h_out = outlet.h;
    (subunits.h*ones(NumberOfSubunits)) = (outlet.h);

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
  end SpeciationOut;
end Processes;
