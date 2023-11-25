within Chemical;
package Processes
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

  model Reaction "Chemical Reaction"
    extends Interfaces.ConditionalKinetics;
    import Chemical.Utilities.Types.InitializationMethods;

    parameter StateSelect n_flowStateSelect = StateSelect.default "State select for n_flow"
      annotation(Dialog(tab="Advanced"));
    parameter InitializationMethods initN_flow = Chemical.Utilities.Types.InitializationMethods.none "Initialization method for n_flow"
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
end Processes;
