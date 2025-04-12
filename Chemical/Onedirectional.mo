within Chemical;
package Onedirectional
  package Processes

    model Reaction "Chemical Reaction"
      extends Interfaces.PartialProcessWithSubstanceData;
      extends Interfaces.ConditionalKinetics(k_forward=1);
      //Real rr_exact2,  kb;
    equation

      rr = kf * Sx * ( 1  -  exp(-duRT));

      /*
    //the same as:
    rr_exact2 = (kf*Sx - kb*Px);
    Kx = kb/kf;
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
</html>",     info="<html>
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

    model ForwardReaction "Chemical Reaction"
      extends Interfaces.PartialProcessWithSubstanceData;
      extends Interfaces.ConditionalKinetics(k_forward=1);

    equation

      rr = kf * Sx;

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
              points={{-60,2},{-60,0},{54,0},{54,0},{18,10},{18,2},{-60,2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-60,-2},{-60,0},{54,0},{54,0},{18,-10},{18,-2},{-60,-2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(revisions="<html>
<p><i>2013-2020 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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
    end ForwardReaction;

    model FastReaction "Chemical Reaction"
      extends Interfaces.PartialProcessWithSubstanceData;

      parameter Real kC = 1;

    equation

      rr = kC * du;



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
</html>",     info="<html>
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
    end FastReaction;

    model FastReactionWithSolutions "Chemical Reaction"
      extends Chemical.Interfaces.PartialProcess;

      parameter stateOfMatter.SubstanceDataParameters productsSubstanceData[nP]
       annotation (choicesAllMatching = true);

      parameter Real kC = 1;

      Interfaces.SolutionPort productSolution[nP] annotation (Placement(transformation(extent={{90,-50},{110,-30}}), iconTransformation(extent={{90,-50},{110,-30}})));
    equation

      rr = kC * du;


      products.definition = productsSubstanceData;

      for i in 1:nP loop
        products[i].solution.T = productSolution[i].T;
        products[i].solution.p = productSolution[i].p;
        products[i].solution.v = productSolution[i].v;
        products[i].solution.n = productSolution[i].n;
        products[i].solution.m = productSolution[i].m;
        products[i].solution.V = productSolution[i].V;
        products[i].solution.G = productSolution[i].G;
        products[i].solution.Q = productSolution[i].Q;
        products[i].solution.I = productSolution[i].I;

        productSolution[i].dH = 0;
        productSolution[i].i = 0;
        productSolution[i].Qj = 0;
        productSolution[i].Ij = 0;
        productSolution[i].nj = 0;
        productSolution[i].mj = 0;
        productSolution[i].Vj = 0;
        productSolution[i].Gj = 0;
        productSolution[i].dV = 0;
      end for;

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
</html>",     info="<html>
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
    end FastReactionWithSolutions;

    model Diffusion "Solute diffusion"
      extends Icons.Diffusion;
      extends Interfaces.PartialChangeSolution(redeclare package stateOfMatterIn = stateOfMatter, redeclare package stateOfMatterOut = stateOfMatter);
      extends Interfaces.ConditionalKinetics(k_forward=1);

      replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
        Interfaces.StateOfMatter "Substance model of inlet"
        annotation (choices(
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

    equation

      n_flow = kf * x_in * ( 1 -  exp(-duRT));


       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end Diffusion;

    model FastDiffusion "Solute diffusion"
      extends Icons.Diffusion;
      extends Interfaces.PartialChangeSolution(redeclare package stateOfMatterIn = stateOfMatter, redeclare package stateOfMatterOut = stateOfMatter);

      replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
        Interfaces.StateOfMatter "Substance model of inlet"
        annotation (choices(
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.Incompressible  "Incompressible"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package stateOfMatterIn =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Real kC = 1;
    equation

      n_flow = kC * du;

       annotation (                 Documentation(revisions="<html>
<p><i>2009-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Diffusion of the substance as equilibration of electro-chemical potentials.</p>
</html>"));
    end FastDiffusion;

    model GasSolubility "Henry's law of gas solubility into liquid."

      extends Icons.GasSolubility;
      extends Interfaces.PartialGasToLiquid;
      extends Interfaces.ConditionalKinetics(k_forward=1);

    equation

      n_flow = kf * x_in * ( 1 -  exp(-duRT));

      annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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

    model GasVolatility "Henry's law of gas volatility from liquid."

      extends Icons.GasSolubility;
      extends Interfaces.PartialLiquidToGas;
      extends Interfaces.ConditionalKinetics(k_forward=1);

    equation

      n_flow = kf * x_in * ( 1 -  exp(-duRT));

      annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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
    end GasVolatility;

    model FastGasSolubility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;
      extends Interfaces.PartialGasToLiquid;


      parameter Real kC = 1;

    equation

      n_flow = kC * du;

      annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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
    end FastGasSolubility;

    model FastGasVolatility "Henry's law of gas solubility in liquid."

      extends Icons.GasSolubility;
      extends Interfaces.PartialLiquidToGas;

      parameter Real kC = 1;

    equation

      n_flow = kC * du;

      annotation (Documentation(revisions="<html>
<p><i>2009-2015 </i></p>
<p><i>by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
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
    end FastGasVolatility;

    model Membrane "Passive transport of the substance through semipermeable membrane"
      extends Icons.Membrane;
      extends Interfaces.PartialChangeSolution;
      extends Interfaces.ConditionalKinetics(k_forward=1);

    equation
      //the main equation

      n_flow = kf * x_in * ( 1 -  exp(-duRT));



      annotation ( Documentation(info="<html>
<p><u><b><font style=\"color: #008000; \">Filtration throught semipermeable membrane.</font></b></u></p>
<p>The penetrating particles are driven by electric and chemical gradient to reach Donnan&apos;s equilibrium.</p>
<p>If zero-flow Donnan&apos;s equilibrium is reached. </p>
</html>", revisions="<html>
<p><i>2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={
            Text(
              extent={{-97,-12},{97,12}},
              textString="%name",
              lineColor={128,0,255},
            origin={69,2},
            rotation=90)}));
    end Membrane;

    model Pump "Prescribed sunstance molar flow"
      extends Interfaces.PartialChangeSolution;
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
              origin={0,-72},
              rotation=360,
              textString="%name")}),        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end Pump;

    model Stream "Flow of whole solution"
      extends Interfaces.PartialChangeSolution;
      extends Boundaries.Internal.ConditionalSolutionFlow;


    equation
      n_flow = c_in * volumeFlow;


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
</html>",     info="<html>
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

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
        Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      /*parameter stateOfMatter.SubstanceDataParameters substanceData
  "Definition of the substance"
    annotation (choicesAllMatching = true);
*/
      parameter stateOfMatter.SubstanceDataParameters subunitsSubstanceData[NumberOfSubunits]
      "Definition of the subunits"
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
      Interfaces.Inlet         inlet annotation (Placement(transformation(extent={{110,-110},{90,-90}}), iconTransformation(extent={{110,-110},{90,-90}})));
      Interfaces.Outlet subunits[NumberOfSubunits] "Subunits of macromolecule" annotation (Placement(transformation(extent={{-56,-14},{-36,66}}),
            iconTransformation(
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
      Chemical.Utilities.Units.URT uRT_out;
      Modelica.Units.SI.ChemicalPotential u_out;
    equation

      inlet.r + inlet.state.u = u_out;

      if NumberOfSubunits>0 then
        (ones(NumberOfSubunits) * subunits.r) = (inlet.r)  -  der(inlet.n_flow)*L;
        for i in 2:NumberOfSubunits loop
          //first subunit is based on inertial potential,
          //other subunits are provided as source
          der(subunits[i].state.u).*TC = subunits[i].r;
        end for;
      end if;


      //amount of macromolecule (all forms in conformation)
      nm*NumberOfSubunits + subunitSolution.nj = 0;

      //change of macromolecule = change of its subunits
      subunits.n_flow = -inlet.n_flow * ones(NumberOfSubunits);

      //mole fraction of all forms in conformation
      xm = nm/solution.n;

      //electrochemical potential of the specific form divided by (Modelica.Constants.R*solution.T)
      uRT_out = log(xm) +
            sum(subunits.state.u/(Modelica.Constants.R*inlet.solution.T) - log(xm)
             * ones(NumberOfSubunits));

      u_out = uRT_out*(Modelica.Constants.R*inlet.solution.T);

      subunits.state.h = (inlet.state.h/NumberOfSubunits)*ones(NumberOfSubunits);

      //inlet.state.u = sum(subunits.state.u - (log(xm)*(Modelica.Constants.R*inlet.solution.T)) * ones(NumberOfSubunits));

      for i in 1:NumberOfSubunits loop
        subunits[i].solution.T = solution.T;
        subunits[i].solution.v = solution.v;
        subunits[i].solution.p = solution.p;
        subunits[i].solution.n = solution.n;
        subunits[i].solution.m = solution.m;
        subunits[i].solution.V = solution.V;
        subunits[i].solution.G = solution.G;
        subunits[i].solution.Q = solution.Q;
        subunits[i].solution.I = solution.I;
      end for;

      //inlet.definition = substanceData;
      subunits.definition = subunitsSubstanceData;

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
</html>",     info="<html>
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

      replaceable package stateOfMatter = Interfaces.Incompressible                    constrainedby
        Interfaces.StateOfMatter
      "Substance model to translate data into substance properties"
         annotation (choicesAllMatching = true);

      parameter stateOfMatter.SubstanceDataParameters substanceData
      "Definition of the substance"
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
      Interfaces.Outlet outletSubstance annotation (Placement(transformation(extent={{90,-110},{110,-90}}), iconTransformation(extent={{90,-110},{110,-90}})));
      Interfaces.Inlet subunits[NumberOfSubunits] "Subunits of macromolecule" annotation (Placement(transformation(extent={{-56,-14},{-36,66}}),
            iconTransformation(
            extent={{10,-40},{-10,40}},
            rotation=90,
            origin={-30,102})));

      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

    //protected

     // stateOfMatter.SubstanceData outletSubstanceDefinition;
      outer DropOfCommons dropOfCommons;
    //  Modelica.Units.SI.MolarEnthalpy h_out;
      Chemical.Utilities.Units.URT uRT_out;
      Modelica.Units.SI.ChemicalPotential u_out; //, u_pure;
    //  Modelica.Units.SI.MolarEntropy s_pure;

    equation

      outletSubstance.state.u = u_out;

      if NumberOfSubunits>0 then
        (ones(NumberOfSubunits) * subunits.r) =(outletSubstance.r) - der(outletSubstance.n_flow)*L;
      end if;

      //amount of macromolecule (all forms in conformation)
      nm*NumberOfSubunits + subunitSolution.nj = 0;

      //change of macromolecule = change of its subunits
      subunits.n_flow =-outletSubstance.n_flow*ones(NumberOfSubunits);

      //mole fraction of all forms in conformation
      xm = nm/solution.n;

      //electrochemical potential of the specific form divided by (Modelica.Constants.R*solution.T)
      uRT_out = log(xm) +
            sum(subunits.state.u/(Modelica.Constants.R*outletSubstance.solution.T) - log(xm)
             * ones(NumberOfSubunits));

      u_out = uRT_out*(Modelica.Constants.R*outletSubstance.solution.T);



      subunits.state.h = (outletSubstance.state.h/NumberOfSubunits)*ones(NumberOfSubunits);

     // h_out =outletSubstance.state.h;
     //  (subunits.state.h*ones(NumberOfSubunits)) =(outletSubstance.state.h);

      outletSubstance.solution.T = solution.T;
      outletSubstance.solution.v = solution.v;
      outletSubstance.solution.p = solution.p;
      outletSubstance.solution.n = solution.n;
      outletSubstance.solution.m = solution.m;
      outletSubstance.solution.V = solution.V;
      outletSubstance.solution.G = solution.G;
      outletSubstance.solution.Q = solution.Q;
      outletSubstance.solution.I = solution.I;

      outletSubstance.definition = substanceData;


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
</html>",     info="<html>
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

  package Boundaries "This package contains boundary models for the stream."
  extends Modelica.Icons.SourcesPackage;

    model Substance "Substance in solution"
      extends Icons.Substance;
      extends Internal.PartialSubstance(m_start=if use_mass_start then mass_start else
         amountOfSubstance_start*stateOfMatter.molarMassOfBaseMolecule(substanceData));

      import Chemical.Utilities.Types.InitializationMethods;


      parameter Boolean use_mass_start=true  "prefere state as mass, otherwise amountOfSubstance"
        annotation (Evaluate=true, choices(checkBox=true));

      parameter Modelica.Units.SI.Mass mass_start=1 "Initial mass of the substance"
        annotation (HideResult=not use_mass_start, Dialog(group="Initialization", enable=use_mass_start));

      parameter Modelica.Units.SI.AmountOfSubstance amountOfSubstance_start=1
      "Initial amount of substance base molecules"
        annotation ( Dialog(group="Initialization", enable=(not use_mass_start)));



      Modelica.Units.SI.Mass mass=n*molarMassOfBaseMolecule "Mass";

      parameter InitializationMethods initAmount = Chemical.Utilities.Types.InitializationMethods.state "Initialization"
        annotation(Dialog(tab="Initialization"));

      parameter Modelica.Units.SI.MolarFlowRate change_start=1
      "Initial molar change of substance base molecules"
        annotation ( Dialog(tab="Initialization", enable=(initAmount == InitializationMethods.derivative)));


    protected
      Modelica.Units.SI.MolarMass molarMassOfBaseMolecule = stateOfMatter.molarMassOfBaseMolecule(substanceDataVar);


      Real logn(stateSelect=StateSelect.prefer, start=log(amountOfSubstance_start))
      "Natural logarithm of the amount of base molecules in solution";

      Real logm(stateSelect=StateSelect.prefer, start=log(mass_start))
      "Natural logarithm of the substance mass in solution";



    initial equation
      if initAmount == InitializationMethods.steadyState then
        r=0;
      elseif initAmount == InitializationMethods.state and not use_mass_start then
        logn=log(amountOfSubstance_start);
      elseif initAmount == InitializationMethods.state and use_mass_start then
        logm=log(mass_start);
      elseif initAmount == InitializationMethods.derivative then
        n_flow = change_start;
      end if;

    equation


      //The main accumulation equation is "der(n)=n_flow"
      // However, the numerical solvers can handle it during equilibration of chemical potential in form of log(n) much better. :-)
      if use_mass_start then
        der(logm) = (n_flow/n) "accumulation of m=exp(logm) [kg]";
      else
        der(logn) = (n_flow/n) "accumulation of n=exp(logn) [mol]";
      end if;
      //der(n) = n_flow;
      n = exp(logn);
      n = exp(logm)*1/molarMassOfBaseMolecule;


        annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-84,22},{92,64}},
              lineColor={128,0,255},
              textString="%name")}), Documentation(revisions="<html>
<p>2009-2025 by Marek Matejak, Ph.D. </p>
</html>",     info="<html>
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

      parameter Interfaces.Incompressible.SubstanceDataParameters substanceData = Chemical.SubstancesOld.Electrone_solid()
                                                                                                                        "Definition of the substance";

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

    model ExternalSubstance "Constant source of molar concentration"
      extends Internal.PartialSubstance(m_start=1, substance(SolutionObserverOnly=true));


      parameter Chemical.Boundaries.Internal.Types.ConcentrationQuantities quantity "Concentration quantity";

      parameter Real FixedValue = 1e-8
      "Fixed value of concentration in selected quantity if useVariableInput=false"
        annotation (HideResult=true, Dialog(enable=not useVariableInput));

      parameter Boolean useVariableInput = false
      "Is amount of substance an input?"
        annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      Modelica.Blocks.Interfaces.RealInput VariableInput=val if useVariableInput
        annotation (HideResult=true, Placement(transformation(extent={{-130,56},{-90,96}})));

      Real value(unit=Chemical.Boundaries.Internal.getUnit(quantity));
    equation

      if not useVariableInput then
        value=FixedValue;
      end if;

      if quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.c_molpm3 then
        value =  1*(substance.x * solutionState.n)/solutionState.V;
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.X_kgpkg then
        value = 1*((substance.x * solutionState.n)/solutionState.m)/stateOfMatter.specificAmountOfParticles(substanceDataVar,
       solutionState.T,
       solutionState.p,
       solutionState.v,
       solutionState.I);
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.b_molpkg then
        value = 1* (substance.x * solutionState.n)/solutionState.m;
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.x_molpmol then
        value = substance.x;
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_Pa then
        value =  1*substance.x*solutionState.p;
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_kPa then
        value*1000 =  substance.x*solutionState.p;
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_mmHg then
        value =  substance.x*solutionState.p * (760/101325);
      elseif quantity == Chemical.Boundaries.Internal.Types.ConcentrationQuantities.p_bar then
        value = 1* Modelica.Units.Conversions.to_bar(substance.x*solutionState.p);
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
    end ExternalSubstance;

    model ExternalIdealGas "Ideal gas substance with defined partial pressure"
      import Chemical;
      extends Chemical.Boundaries.Internal.PartialSubstanceOld
                                                           (
        redeclare package stateOfMatter = gasModel "Ideal Gas from MSL",
        m_start=1,
        substance(SolutionObserverOnly=true));

       replaceable package gasModel = Chemical.Interfaces.IdealGasMSL constrainedby
        Chemical.Interfaces.StateOfMatter "Gas substance model"
        annotation (choices(
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGas        "Ideal Gas"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasMSL     "Ideal Gas from MSL"),
          choice(redeclare package gasModel =
            Chemical.Interfaces.IdealGasShomate "Ideal Gas using Shomate model")));

      parameter Boolean usePartialPressureInput = false
      "=true, if fixed partial pressure is from input instead of parameter"
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

      parameter Modelica.Units.SI.Pressure PartialPressure=1e-05
        "Fixed partial pressure if usePartialPressureInput=false" annotation (
         HideResult=true, Dialog(enable=not usePartialPressureInput));

      Modelica.Blocks.Interfaces.RealInput partialPressure(start=
            PartialPressure, final unit="Pa")=p if usePartialPressureInput
      "Partial pressure of gas = total pressure * gas fraction"
        annotation (HideResult=true,Placement(transformation(extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,72}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,72})));

      Modelica.Units.SI.Pressure p "Current partial pressure";

    equation

      if not usePartialPressureInput then
        p=PartialPressure;
      end if;

      //mole fraction
      substance.x = p / solutionState.p;



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

    model TerminalInflow "Molar pump of substance to system"
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

      outlet.n_flow = -q;

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
    end TerminalInflow;

    model TerminalOutflow "Molar pump of substance out of system"
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
    end TerminalOutflow;

    model Clearance "ClearanceFlow of whole solution"
      extends Internal.PartialSubstanceInlet;

      parameter Modelica.Units.SI.VolumeFlowRate Clearance
      "Physiological clearance of the substance";

    equation

      assert(Clearance>=-Modelica.Constants.eps, "Clearance can not be negative!");

      inlet.n_flow = substance.c * Clearance;

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
</html>",     info="<html>
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
      extends Internal.PartialSubstanceInlet;

      parameter Modelica.Units.SI.Time HalfTime
      "Degradation half time. The time after which will remain half of initial concentration in the defined volume when no other generation, clearence and degradation exist.";

      Modelica.Units.SI.AmountOfSubstance n;

    equation

      n = substance.x*inlet.solution.n;
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
</html>",     info="<html>
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

    package Internal "Partials and Internal functions"
    extends Modelica.Icons.InternalPackage;

      partial model PartialBoundary
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

        parameter stateOfMatter.SubstanceDataParameters substanceData    "Definition of the substance"
          annotation (choicesAllMatching = true, Dialog(enable = not useInlet));

        parameter Boolean useInlet=false   "Use inlet conector?"
            annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L
         annotation( Dialog(tab = "Advanced"));

        parameter Modelica.Units.SI.MolarFlowRate n_flow_assert(max=0) = -dropOfCommons.n_flow_reg "Assertion threshold for negative molar flows"
          annotation(Dialog(tab="Advanced"));


       outer Modelica.Fluid.System system "System wide properties";

        Chemical.Interfaces.Inlet inlet(
          redeclare package stateOfMatter=stateOfMatter,
          r=r_in,
          n_flow=n_flow_in,
          solution=solutionState) if useInlet
          "The substance entering"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                                                                            iconTransformation(extent={{-110,-10},{-90,10}})));

        Chemical.Interfaces.Outlet outlet(
          redeclare package stateOfMatter = stateOfMatter,
          r=r_out,
          n_flow=n_flow_out,
          state=state_out,
          solution=solutionState,
          definition=substanceDataVar)
                     "The substance exiting" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

       stateOfMatter.SubstanceData substanceDataVar;

       Modelica.Units.SI.MolarFlowRate n_flow "Total molar change of substance";

       Modelica.Units.SI.EnthalpyFlowRate h_flow "Enthalpy change";

       Modelica.Units.SI.AmountOfSubstance n
          "Amount of base molecules inside all clusters in compartment";

      protected

        parameter Modelica.Units.SI.Mass m_start "Start value for mass of the substance";

        outer Chemical.DropOfCommons dropOfCommons;

        stateOfMatter.InputSubstanceData substanceDefinition;
        Chemical.Interfaces.SolutionState solutionState;

        Modelica.Units.SI.ChemicalPotential r = state_out.u - state_in.u;

        Modelica.Units.SI.ChemicalPotential r_in;
        Modelica.Units.SI.ChemicalPotential r_out;
        Modelica.Units.SI.MolarFlowRate n_flow_in;
        Modelica.Units.SI.MolarFlowRate n_flow_out;

        Chemical.Interfaces.InputSubstanceState state_in;
        Chemical.Interfaces.SubstanceState state_out;

      equation

       der(n_flow_out)*L = r_out;
       der(n_flow_in)*L = r_in - r;

       connect(substanceDefinition,inlet.definition);
        substanceDataVar = substanceDefinition;
        connect(state_in,inlet.state);

        if not useInlet then
          r_in=0;
          n_flow_in=0;

          state_in.h=0;
          substanceDataVar = substanceData;
        end if;

        n_flow = n_flow_in + n_flow_out;
        h_flow = n_flow_out*state_out.h + n_flow_in*state_in.h;

       annotation (
         Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Ph.D.</p>
</html>"));
      end PartialBoundary;

      partial model PartialSubstance
        import Chemical;

        extends PartialBoundary;

        parameter Boolean useSolution = true "Use solution connector?"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Conditional inputs"));

        parameter Chemical.Interfaces.SolutionStateParameters solutionParam "Constant chemical solution state if not from rear or input"
          annotation (Dialog(enable=not useSolution and not useRear));

        Chemical.Interfaces.SolutionPort solution(
            T=solutionPortState.T,
            p=solutionPortState.p,
            v=solutionPortState.v,
            n=solutionPortState.n,
            m=solutionPortState.m,
            V=solutionPortState.V,
            G=solutionPortState.G,
            Q=solutionPortState.Q,
            I=solutionPortState.I,
            i=substance.i,
            dH=substance.dH,
            dV=substance.dV,
            nj=substance.nj,
            mj=substance.mj,
            Vj=substance.Vj,
            Gj=substance.Gj,
            Qj=substance.Qj,
            Ij=substance.Ij)
              if useSolution "To connect substance with solution, where is pressented"
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}), iconTransformation(extent={{-70,-110},{-50,-90}})));


       stateOfMatter.BaseProperties substance(
         substanceDataVar=substanceDataVar,
         solutionState=solutionState,
         FixedSubstanceData=not useInlet,
         substanceData=substanceData,
         amountOfBaseMolecules=n,
         m_start=m_start,
         n_flow=n_flow,
         h_flow=h_flow);


      protected

           Chemical.Interfaces.SolutionState solutionPortState;

      equation

       //electro-chemical potential of the substance in the solution
       state_out.u = substance.u;
       state_out.h = substance.h;

       if (useSolution and not useInlet) or (not useSolution) then
          solutionState=solutionPortState;
        end if;

        if not useSolution and not useInlet then
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



       annotation (
         Documentation(revisions="<html>
<p><i>2009-2025</i></p>
<p>Marek Matejak, Ph.D.</p>
</html>"));
      end PartialSubstance;

      partial model PartialSubstanceInlet

       outer Modelica.Fluid.System system "System wide properties";

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

        Interfaces.Inlet inlet(redeclare package stateOfMatter = stateOfMatter) "The substance"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));

        stateOfMatter.BaseProperties substance(
          SolutionObserverOnly=true,
          substanceDataVar=inlet.definition,
          solutionState=inlet.solution,
          FixedSubstanceData=false,
          m_start=1,
          n_flow=0,
          h_flow=0);

      equation

       substance.u = inlet.state.u;

       annotation (
         Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end PartialSubstanceInlet;

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

</html>",   info="<html>
<p>The boundaries are Sorces and Sinks, as well as Volumes, that are conceptually a source and a sink with extra equations and act as loop breakers in closes cycles, and therefore are also boundaries.</p>
</html>"));
  end Boundaries;

  package Sensors "Chemical sensors"
    extends Modelica.Icons.SensorsPackage;

    model MolarFlowSensor "Measure of molar flow"
      extends Modelica.Icons.RoundSensor;

      Modelica.Blocks.Interfaces.RealOutput molarFlowRate(final unit="mol/s") annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-100})));

      Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
      Interfaces.Outlet outlet annotation (Placement(transformation(extent={{92,-10},{112,10}})));
    equation
      molarFlowRate = inlet.n_flow;

      connect(inlet, outlet) annotation (Line(
          points={{-98,0},{102,0}},
          color={158,66,200},
          thickness=0.5));
     annotation (
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={
            Line(
              points={{70,-10},{90,-10}},
              color={127,0,127}),
            Line(
              points={{70,10},{90,10}},
              color={127,0,127}),
            Line(
              points={{-90,10},{-70,10}},
              color={127,0,127}),
            Line(
              points={{-90,-10},{-70,-10}},
              color={127,0,127}),
            Text(
              extent={{-31,-5},{28,-64}},
              lineColor={0,0,0},
              textString="dn")}));
    end MolarFlowSensor;

    model MoleFractionSensor "Measure of mole fraction"
      extends Modelica.Icons.RoundSensor;
      extends Internal.SubstanceSensor;

      Modelica.Blocks.Interfaces.RealOutput moleFraction(final unit="1")
      "Mole fraction of the substance"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={100,0},
          rotation=0)));

    equation

      moleFraction = x;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
              textString="x"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MoleFractionSensor;

    model MolalitySensor "Measure of molality of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Internal.SubstanceSensor;

       Modelica.Blocks.Interfaces.RealOutput molality(final unit="mol/kg")
      "Molality of the substance (amount of substance per mass of solvent)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={100,0},
          rotation=0)));

    equation

      molality = x * (inlet.solution.n/inlet.solution.m);

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="b"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolalitySensor;

    model MolarConcentrationSensor "Measure of molarity of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Internal.SubstanceSensor;


       Modelica.Blocks.Interfaces.RealOutput molarConcentration(final unit="mol/m3", displayUnit="mol/l")
      "Molarity of the substance (amount of substance in one liter of whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={100,0},
          rotation=0)));

    equation

      molarConcentration = x * (inlet.solution.n/inlet.solution.V);

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="c"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MolarConcentrationSensor;

    model MassFractionSensor "Measure of mass fraction of the substance"
      extends Modelica.Icons.RoundSensor;
      extends Internal.SubstanceSensor;


       Modelica.Blocks.Interfaces.RealOutput massFraction(final unit="kg/kg")
      "Mass fraction of the substance (mass of the substance per mass of the whole solution)"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={100,0},
          rotation=0)));

    equation

      massFraction = (x / stateOfMatter.specificAmountOfParticles(substanceData)) * (inlet.solution.n/inlet.solution.m);

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="mx"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end MassFractionSensor;

    model PartialPressureSensor
    "Measure of partial pressure of the substance in gaseous solution"
      extends Modelica.Icons.RoundSensor;
      extends Internal.SubstanceSensor;

       Modelica.Blocks.Interfaces.RealOutput partialPressure(final unit="Pa")
      "Partial pressure of the substance in gaseous solution"
       annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,-60}), iconTransformation(
            extent={{-20,-20},{20,20}},
            origin={100,0},
          rotation=0)));

    equation

      partialPressure = x*inlet.solution.p;

     annotation (
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}),   graphics={
            Text(
              extent={{-31,-3},{28,-62}},
              lineColor={0,0,0},
            textString="p"),
            Line(
              points={{70,0},{80,0}},
              color={127,0,127})}),
        Documentation(revisions="<html>
<p><i>2009-2015</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end PartialPressureSensor;

    package Internal

      model SubstanceSensor "Base class for sensor based on inlet substance and solution properties"
        extends Boundaries.Internal.PartialSubstanceInlet;
      equation
        inlet.n_flow=0;
      end SubstanceSensor;

      package Types
        type Quantities = enumeration(
            c_molpm3 "Concentration (mmol/L)",
            X_kgpkg "Mass fraction (kg/kg)",
            b_molpkg "Molality (mol/kg)",
            x_molpmol "Mole fraction (mol/mol)",
            p_Pa "Partial pressure (Pa)",
            p_mmHg "Partial pressure (mmHg)",
            p_bar "Partial pressure (bar)",
            u_Jpmol "Steadystate electro-chemical potential (J/mol)",
            u_kJpmol "Steadystate electro-chemical potential (kJ/mol)",
            r_Jpmol "Inertial electro-chemical potential (J/mol)",
            r_kJpmol "Inertial electro-chemical potential (kJ/mol)",
            u_total_Jpmol "Total electro-chemical potential (J/mol)",
            u_total_kJpmol "Total electro-chemical potential (kJ/mol)",
            h_Jpmol "Specific enthalpy (J/mol)",
            s_JpmolK "Specific enthropy (J/(mol.K))");
        type InitializationModelSensor = enumeration(
          steadyState
            "Steady state initialization (derivatives of states are zero)",
          state
            "Initialization with initial output state") "Initialization modes for sensor lowpass";
      end Types;

      function getQuantity "Computes selected quantity from state"
        extends Modelica.Icons.Function;

        replaceable package stateOfMatter = Interfaces.Incompressible constrainedby
          Interfaces.StateOfMatter
        "Substance model to translate data into substance properties"
          annotation (choicesAllMatching=true,
            Documentation(info="<html>
      <p>Medium Model for the function. Make sure it implements the needed functions.</p>
        </html>"));

        input Modelica.Units.SI.ChemicalPotential u "Electro-chemical potential";
        input Modelica.Units.SI.MolarEnthalpy h "Molar enthalpy";
        input Modelica.Units.SI.ChemicalPotential r "Inertial electro-chemical potential";
        input Types.Quantities quantity "What to measure?";
        input stateOfMatter.SubstanceData substanceData "Data record of substance";
        input Modelica.Units.SI.Temperature temperature=298.15 "Temperature";
        input Modelica.Units.SI.Pressure pressure=100000 "Pressure";
        input Modelica.Units.SI.ElectricPotential electricPotential=0
         "Electric potential of the substance";
        input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0
         "Ionic strengh (mole fraction based)";
        input Modelica.Units.SI.Mass solutionMass;
        input Modelica.Units.SI.AmountOfSubstance solutionAmount;
        input Modelica.Units.SI.Volume solutionVolume;

        output Real value;

        /*
    c_molpm3 "Concentration (mmol/L)",
    X_kgpkg "Mass fraction (kg/kg)",
    b_molpkg "Molality (mol/kg)",
    x_molpmol "Mole fraction (mol/mol)",
    p_Pa "Partial pressure (Pa)",
    p_mmHg "Partial pressure (mmHg)",
    p_bar "Partial pressure (bar)",
    u_Jpmol "Steadystate electro-chemical potential (J/mol)",
    u_kJpmol "Steadystate electro-chemical potential (kJ/mol)",
    r_Jpmol "Inertial electro-chemical potential (J/mol)",
    r_kJpmol "Inertial electro-chemical potential (kJ/mol)",
    u_total_Jpmol "Total electro-chemical potential (J/mol)",
    u_total_kJpmol "Total electro-chemical potential (kJ/mol)",
    h_Jpmol "Specific enthalpy (J/mol)",
    s_JpmolK "Specific enthropy (J/(mol.K))"
    */

      protected
        Modelica.Units.SI.ChargeNumberOfIon z;
        Modelica.Units.SI.ChemicalPotential u_Pure;
        Modelica.Units.SI.MoleFraction a,x;
        Modelica.Units.SI.ActivityCoefficient gamma
        "Activity coefficient of the substance";

      algorithm
        z :=  stateOfMatter.chargeNumberOfIon( substanceData,
         temperature,
         pressure,
         electricPotential,
         moleFractionBasedIonicStrength);
        gamma := stateOfMatter.activityCoefficient(substanceData,temperature,pressure,electricPotential,moleFractionBasedIonicStrength);


        u_Pure := stateOfMatter.chemicalPotentialPure(
         substanceData,
         temperature,
         pressure,
         electricPotential,
         moleFractionBasedIonicStrength)
         + z*Modelica.Constants.F*electricPotential;

        a := exp((u - u_Pure)/(Modelica.Constants.R*temperature));
        x := a/gamma;

        if quantity == Types.Quantities.c_molpm3 then
          value := (x * solutionAmount)/solutionVolume;
        elseif quantity == Types.Quantities.X_kgpkg then
          value := ((x * solutionAmount)/solutionMass)/stateOfMatter.specificAmountOfParticles(substanceData,
         temperature,
         pressure,
         electricPotential,
         moleFractionBasedIonicStrength);
        elseif quantity == Types.Quantities.b_molpkg then
          value := (x * solutionAmount)/solutionMass;
        elseif quantity == Types.Quantities.x_molpmol then
          value := x;
        elseif quantity == Types.Quantities.p_Pa then
          value := x*pressure;
        elseif quantity == Types.Quantities.p_mmHg then
          value := x*pressure * (760/101325);
        elseif quantity == Types.Quantities.p_bar then
          value := Modelica.Units.Conversions.to_bar(x*pressure);
        elseif quantity == Types.Quantities.u_Jpmol then
          value := u;
        elseif quantity == Types.Quantities.u_kJpmol then
          value := u/1000;
        elseif quantity == Types.Quantities.h_Jpmol then
          value := h;
        elseif quantity == Types.Quantities.s_JpmolK then
          value := (h-u)/temperature;
        else
          value :=0;
        end if;

        annotation (Documentation(info="<html>
<p>Helper function to get a quantity from an Thermofluid state.</p>
</html>"));
      end getQuantity;

      function getUnit "Returns unit of input quantity"
        extends Modelica.Icons.Function;

        input Types.Quantities quantity;
        output String unit;

      algorithm

        if quantity == Types.Quantities.c_molpm3 then
          unit := "mol/m3";
        elseif quantity == Types.Quantities.X_kgpkg then
          unit := "kg/kg";
        elseif quantity == Types.Quantities.b_molpkg then
          unit := "mol/kg";
        elseif quantity == Types.Quantities.x_molpmol then
          unit := "mol/mol";
        elseif quantity == Types.Quantities.p_Pa then
          unit := "Pa";
        elseif quantity == Types.Quantities.p_mmHg then
          unit := "mmHg";
        elseif quantity == Types.Quantities.p_bar then
          unit := "bar";
        elseif quantity == Types.Quantities.u_Jpmol then
          unit := "J/mol";
        elseif quantity == Types.Quantities.u_kJpmol then
          unit := "kJ/mol";
        elseif quantity == Types.Quantities.h_Jpmol then
          unit :="J/mol";
        elseif quantity == Types.Quantities.s_JpmolK then
          unit := "J/(mol.K)";
        else
          unit :="";
        end if;

        annotation (Documentation(info="<html>
<p>Helper function to get the unit for a quantity.</p>
</html>"));
      end getUnit;
    end Internal;

    model SingleSensorSelect "Sensor with selectable measured quantity"
      import Chemical.Sensors.Internal.Types.Quantities;
      import InitMode = Chemical.Sensors.Internal.Types.InitializationModelSensor;
      extends Internal.SubstanceSensor;

      parameter Integer digits(min=0) = 1 "Number of displayed digits";
      parameter Quantities quantity "Quantity the sensor measures";
      parameter Boolean outputValue = false "Enable sensor-value output"
        annotation(Dialog(group="Output Value"));
      parameter Boolean filter_output = false "Filter sensor-value to break algebraic loops"
        annotation(Dialog(group="Output Value", enable=outputValue));
      parameter InitMode init=InitMode.steadyState "Initialization mode for sensor lowpass"
        annotation(Dialog(tab="Initialization", enable=filter_output));
      parameter Real value_0(unit=Internal.getUnit(quantity)) = 0 "Initial output state of sensor"
        annotation(Dialog(tab="Initialization", enable=filter_output and init == InitMode.state));
      parameter Modelica.Units.SI.Time TC=0.1 "PT1 time constant" annotation (Dialog(tab="Advanced", enable=outputValue and filter_output));

      Modelica.Blocks.Interfaces.RealOutput value_out(unit=Internal.getUnit(quantity)) = value if outputValue "Measured value [variable]"
        annotation (Placement(transformation(extent={{80,-20},{120,20}})));

      output Real value(unit=Internal.getUnit(quantity)) "Computed value of the selected quantity";

    protected
      outer DropOfCommons dropOfCommons;

      Real direct_value(unit=Internal.getUnit(quantity));

      function getQuantity = Internal.getQuantity (redeclare package stateOfMatter =
              stateOfMatter)                                                       "Quantity compute function"
        annotation (Documentation(info="<html>
      <p>This function computes the selected quantity from state. r and rho_min are neddet for the quantities r/p_total and v respectively.</p>
      </html>"));

    initial equation
      if filter_output and init==InitMode.steadyState then
        value= direct_value;
      elseif filter_output then
        value = value_0;
      end if;

    equation

      direct_value = getQuantity(inlet.state.u, inlet.state.h, inlet.r, quantity, substanceData, temperature, pressure, electricPotential, moleFractionBasedIonicStrength,
       inlet.solution.m, inlet.solution.n, inlet.solution.V);


      if filter_output then
        der(value) * TC = direct_value-value;
      else
        value = direct_value;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-54,24},{66,-36}},
              lineColor={0,0,0},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-100,0},{0,0}},
              color={28,108,200},
              thickness=0.5),
            Rectangle(
              extent={{-60,30},{60,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,30},{60,-30}},
              textColor={28,108,200},
              textString=DynamicSelect("value", String(
                  value,
                  format="1."+String(digits)+"f"))),
            Text(
              extent={{0,25},{60,75}},
              textColor={175,175,175},
              textString="%quantity")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>Sensor for measuring a selectable quantity.</p>
<p>This sensor can be connected to a fluid stream without a junction.</p>
</html>"));
    end SingleSensorSelect;
  end Sensors;

  package Topology
    extends Modelica.Icons.Package;

    model SplitterT1 "Splits a flow into two subflows"

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

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Interfaces.Inlet inlet(redeclare package stateOfMatter = stateOfMatter) annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,0})));
      Interfaces.Outlet outletA(redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,100})));
      Chemical.Interfaces.Outlet outletB(redeclare package stateOfMatter = stateOfMatter)
        annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,-100})));
      SplitterN splitterN(redeclare package stateOfMatter = stateOfMatter,final N=2, final L=L)
        annotation (Placement(transformation(extent={{-28,-10},{-8,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(splitterN.inlet, inlet) annotation (Line(
          points={{-28,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet[2], outletA) annotation (Line(
          points={{-8,0.5},{-8,50},{0,50},{0,100}},
          color={158,66,200},
          thickness=0.5));
      connect(outletB, splitterN.outlet[1]) annotation (Line(
          points={{0,-100},{0,-0.5},{-8,-0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,100},{60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{20,-60},{60,-100}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SplitterT1;

    model SplitterT2 "Splits a flow into two subflows"

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Interfaces.Inlet inletProcess annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
      Interfaces.Outlet outletSubstance annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,100})));
      Interfaces.Outlet outletSubstance1 annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));
      SplitterN splitterN(final N=2, final L=L)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(splitterN.inlet, inletProcess) annotation (Line(
          points={{-40,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet[1], outletSubstance1) annotation (Line(
          points={{-20,-0.5},{-20,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(outletSubstance, splitterN.outlet[2]) annotation (Line(
          points={{0,100},{0,0.5},{-20,0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{60,-20},{100,-60}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SplitterT2;

    model JunctionT1 "2 to 1 T-Junction"

      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,100})));
      Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,-100})));
      JunctionN junctionN(final N=2, final L=L,
         final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-34,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionN.outlet, outlet) annotation (Line(
          points={{-44,1.33227e-15},{-72,1.33227e-15},{-72,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(inletA, junctionN.inlets[1]) annotation (Line(
          points={{0,100},{0,0.5},{-24,0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionN.inlets[2], inletB) annotation (Line(
          points={{-24,-0.5},{0,-0.5},{0,-100}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{20,100},{60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{20,-60},{60,-100}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionT1;

    model JunctionT2 "2 to 1 T-Junction"

      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,100})));
      Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={100,0})));
      JunctionN junctionN(final N=2, final L=L, final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-20,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionN.inlets[2], inletB) annotation (Line(
          points={{-10,-0.5},{36,-0.5},{36,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(inletA, junctionN.inlets[1]) annotation (Line(
          points={{0,100},{0,0.5},{-10,0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionN.outlet, outlet) annotation (Line(
          points={{-30,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{80,-20},{120,-60}},
              textColor={175,175,175},
              textString="B")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionT2;

    model SplitterX "Splits a flow into three subflows"

      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Inlet inlet annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
      Chemical.Interfaces.Outlet outletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,100})));
      Chemical.Interfaces.Outlet outletB
        annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={3.55271e-15,-100})));
      Chemical.Interfaces.Outlet outletC annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));
      SplitterN splitterN(final N=3, final L=L)
        annotation (Placement(transformation(extent={{-32,-10},{-12,10}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(splitterN.inlet, inlet) annotation (Line(
          points={{-32,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet[3], outletA) annotation (Line(
          points={{-12,0.666667},{0,0.666667},{0,100}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet[2], outletC) annotation (Line(
          points={{-12,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet[1], outletB)
        annotation (Line(
          points={{-12,-0.666667},{0,-0.666667},{0,-100},{3.55271e-15,-100}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-20,100},{-60,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{20,-60},{60,-100}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{50, 60},{90, 20}},
              textColor={175,175,175},
              textString="C")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SplitterX;

    model JunctionX1 "2 to 2 X-Junction"

      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Outlet outleta annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Outlet outletb annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,-100})));
      Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,100})));
      Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={100,0})));
      JunctionNM junctionNM(N=2, M=2, final L=L, final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation
      connect(inletB, junctionNM.inlets[2]) annotation (Line(
          points={{100,0},{54,0},{54,-0.5},{10,-0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionNM.outlets[1], outleta) annotation (Line(
          points={{-10,0.5},{-54,0.5},{-54,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(inletA, junctionNM.inlets[1]) annotation (Line(
          points={{0,100},{0,40},{40,40},{40,0.5},{10,0.5}},
          color={158,66,200},
          thickness=0.5));
      connect(outletb, junctionNM.outlets[2]) annotation (Line(
          points={{0,-100},{0,-40},{-40,-40},{-40,-0.5},{-10,-0.5}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{-60,-20},{-100,-60}},
              textColor={175,175,175},
              textString="a"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="b")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>Medium package used in the component. Make sure it is the same one as all the components connected to all fluid ports are using. </p>
</html>"));
    end JunctionX1;

    model JunctionX2 "2 to 2 X-Junction"

      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Outlet outleta annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Outlet outletb annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));
      Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,100})));
      Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,-100})));
      JunctionNM junctionNM(N=2, M=2, final L=L, final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,20})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionNM.outlets[1], outleta) annotation (Line(
          points={{-0.5,10},{-1,10},{-1,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionNM.outlets[2], outletb) annotation (Line(
          points={{0.5,10},{1,10},{1,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionNM.inlets[1], inletA) annotation (Line(
          points={{-0.5,30},{-0.5,66},{0,66},{0,100}},
          color={158,66,200},
          thickness=0.5));
      connect(inletB, junctionNM.inlets[2]) annotation (Line(
          points={{0,-100},{0,-40},{40,-40},{40,52},{0.5,52},{0.5,30}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,-100},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{-60,-20},{-100,-60}},
              textColor={175,175,175},
              textString="a"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="b")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionX2;

    model JunctionX3 "3 to 1 X-Junction"

      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Outlet outlet annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.Inlet inletA annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,100})));
      Chemical.Interfaces.Inlet inletB annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=180,
            origin={100,0})));
      Chemical.Interfaces.Inlet inletC annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,-100})));
      JunctionN junctionN(final N=3, final L=L,
        final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionN.inlets[2], inletB) annotation (Line(
          points={{-30,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(inletC, junctionN.inlets[3]) annotation (Line(
          points={{0,-100},{0,-0.666667},{-30,-0.666667}},
          color={158,66,200},
          thickness=0.5));
      connect(inletA, junctionN.inlets[1]) annotation (Line(
          points={{0,100},{0,0.666667},{-30,0.666667}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionN.outlet, outlet) annotation (Line(
          points={{-50,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{0,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,100}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{0,-100}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-60,100},{-20,60}},
              textColor={175,175,175},
              textString="A"),
            Text(
              extent={{50,20},{90,60}},
              textColor={175,175,175},
              textString="B"),
            Text(
              extent={{60,-100},{20,-60}},
              textColor={175,175,175},
              textString="C")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionX3;

    model SplitterN "Splitter with one inlet and N outlets"

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

      parameter Integer N(min=1) = 1 "Number of outputs";
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Interfaces.Inlet inlet(redeclare package stateOfMatter = stateOfMatter) "inlet" annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));
      Interfaces.Outlet outlet[N](redeclare package stateOfMatter = stateOfMatter) "vector of N outlets"
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential r_mix;

    equation
      der(inlet.n_flow)*L = inlet.r - r_mix;

      for i in 1:N loop
        der(outlet[i].n_flow)*L = outlet[i].r - r_mix;
        outlet[i].definition = inlet.definition;
        outlet[i].solution = inlet.solution;

        outlet[i].state.u = inlet.state.u;
        outlet[i].state.h = inlet.state.h;
      end for;

      sum(outlet.n_flow) + inlet.n_flow = 0;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{0,0},{96,10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{96,-10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{-100,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{90,80},{50,40}},
              textColor={175,175,175},
              textString="%N")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SplitterN;

    model JunctionN "Junction with N inlets and one outlet"

      parameter Integer N(min=1) = 1 "Number of inlets";
      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small molar flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Inlet inlets[N] "vector of N inlets"
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,-20},{-80,20}})));
      Chemical.Interfaces.Outlet outlet "outlet"
        annotation (Placement(transformation(extent={{80,-20},{120,20}}), iconTransformation(extent={{80,-20},{120,20}})));

      // these are needed by DynamicJunctionN
      output Real w[N](each unit="1") "regularized weighting factor for specific enthalpy";

    protected
      outer Chemical.DropOfCommons dropOfCommons;

      Modelica.Units.SI.ChemicalPotential u[N]=inlets.state.u "(steady molar-flow) electro-chemical potential";
      Modelica.Units.SI.MolarEnthalpy h[N]=inlets.state.h "molar enthapy at inlets";

      Modelica.Units.SI.ChemicalPotential u_mix "(steady mass-flow) electro-chemical potential at the outlet";

      Modelica.Units.SI.ChemicalPotential r_mix "inertial electro-chemical potential divided by R*T at outlet";
      Modelica.Units.SI.MolarEnthalpy h_mix "molar enthalpy at outlet";

      Modelica.Units.SI.ChemicalPotential r_in[N];

    equation
      sum(inlets.n_flow) + outlet.n_flow = 0;

      for i in 1:N loop
        der(inlets[i].n_flow) * L = inlets[i].r - r_in[i];

        u[i] + r_in[i] = u_mix + r_mix;
        w[i] = (abs(inlets[i].n_flow)+n_flow_eps) / (sum(abs(inlets.n_flow))+N*n_flow_eps);

        //assert(outlet.definition.DGH >= inlets[i].definition,"Topology connections allowed only for the same substance.");
        //assert(outlet.solution.T == inlets[i].solution.T,"Topology connections allowed only for the same chemical solution.");
      end for;
      der(outlet.n_flow) * L =  outlet.r - r_mix;

      u_mix = sum(w.*u);
      h_mix = sum(w.*h);

      outlet.state.u = u_mix;
      outlet.state.h = h_mix;

      outlet.definition = inlets[1].definition;
      outlet.solution = inlets[1].solution;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{-100,10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{-100,-10}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-90,80},{-50,40}},
              textColor={175,175,175},
              textString="%N")}),
          Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionN;

    model JunctionNM "Junction with N inlets and M outlets"

      parameter Integer N(min=1) = 1 "Number of inputs";
      parameter Integer M(min=1) = 1 "Number of outputs";
      parameter Modelica.Units.SI.MolarFlowRate n_flow_eps=dropOfCommons.n_flow_reg "Regularization threshold for small mass flows"
        annotation (Dialog(tab="Advanced"));
      parameter Chemical.Utilities.Units.Inertance L=dropOfCommons.L "Inertance on each Branch of Component" annotation (Dialog(tab="Advanced"));

      Chemical.Interfaces.Inlet inlets[N] "vector of N inlets" annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,0}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-100,0})));
      Chemical.Interfaces.Outlet outlets[M] "vector of N outlets" annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={100,0})));
      SplitterN splitterN(
        final N=M,
        final L=L)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={18,0})));
      JunctionN junctionN(
        final N=N,
        final L=L,
        final n_flow_eps=n_flow_eps)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-14,0})));

    protected
      outer Chemical.DropOfCommons dropOfCommons;

    equation

      connect(junctionN.inlets, inlets) annotation (Line(
          points={{-24,0},{-62,0},{-62,0},{-100,0}},
          color={158,66,200},
          thickness=0.5));
      connect(junctionN.outlet, splitterN.inlet) annotation (Line(
          points={{-4,0},{8,0}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterN.outlet, outlets) annotation (Line(
          points={{28,0},{100,0}},
          color={158,66,200},
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Line(
              points={{-100,0},{100,0}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{-100,10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{-100,-10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{96,10}},
              color={158,66,200},
              thickness=0.5),
            Line(
              points={{0,0},{96,-10}},
              color={158,66,200},
              thickness=0.5),
            Ellipse(
              extent={{-6,6},{6,-6}},
              lineColor={158,66,200},
              fillColor={194,138,221},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Text(
              extent={{-90,80},{-50,40}},
              textColor={175,175,175},
              textString="%N"),
            Text(
              extent={{90,80},{50,40}},
              textColor={175,175,175},
              textString="%M")}),
          Diagram(coordinateSystem(preserveAspectRatio=false)));
    end JunctionNM;
  annotation (Documentation(revisions="<html>
<p><img src=\"modelica:/Chemical/Resources/dlr_logo.png\"/>(c) 2020-2021, DLR, Institute of System Dynamics and Control</p>

</html>"),   Icon(graphics={
          Line(
            points={{-80,0},{12,0}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{12,0},{80,-80}},
            color={158,66,200},
            thickness=0.5),
          Line(
            points={{12,0},{80,80}},
            color={158,66,200},
            thickness=0.5),
          Ellipse(
            extent={{6,6},{18,-6}},
            lineColor={158,66,200},
            fillColor={194,138,221},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}));
  end Topology;

  package Examples "Examples that demonstrate usage of chemical library"
  extends Modelica.Icons.ExamplesPackage;

    model SimpleReaction
      "The simple chemical reaction A<->B with equilibrium B/A = 2"
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real K = 2 "Dissociation constant of the reaction";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Onedirectional.Boundaries.Substance
                                    A(use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-52,-8},{-32,12}})));

      Chemical.Onedirectional.Processes.Reaction
                                  reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(K))},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
      Chemical.Onedirectional.Boundaries.Substance
                                    B(
        useInlet=true,
        use_mass_start=false,
        amountOfSubstance_start=0.1)
                       annotation (Placement(transformation(extent={{42,-8},{62,12}})));

      inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{58,64},{78,84}})));
    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-48,-8},{-48,-92},{60,-92},{60,-98}},
          color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{46,-8},{46,-92},{60,-92},{60,-98}},
                                         color={127,127,0}));
      connect(A.outlet, reaction2_1.substrates[1]) annotation (Line(
          points={{-32,2},{-10,2}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], B.inlet) annotation (Line(
          points={{10,2},{42,2}},
          color={158,66,200},
          thickness=0.5));
      annotation (Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A and substance B, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that mole fraction (A.x and B.x) are always summed to 1 for the solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end SimpleReaction;

    model SimpleReaction2
      "The simple chemical reaction A+B<->C with equilibrium [C]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
      import Chemical;
       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-32,4},{-12,24}})));
      Chemical.Processes.Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=2, DfG=-R*T_25degC*log(Kx))},
        nS=2,
        nP=1) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Boundaries.Substance B(
        substanceData(MolarWeight=1),
        use_mass_start=false,
        amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Boundaries.Substance C(
        use_mass_start=false,
        amountOfSubstance_start=0.1,
        useInlet=true) annotation (Placement(transformation(extent={{50,-8},{70,12}})));

    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-28,4},{-28,-90},{60,-90},{60,-98}},
          color={127,127,0}));
      connect(C.solution, solution.solution) annotation (Line(points={{54,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                 color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
            {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

      connect(B.outlet, reaction2_1.substrates[1]) annotation (Line(
          points={{-14,-14},{-4,-14},{-4,1.75},{4,1.75}},
          color={158,66,200},
          thickness=0.5));
      connect(A.outlet, reaction2_1.substrates[2]) annotation (Line(
          points={{-12,14},{-4,14},{-4,2.25},{4,2.25}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], C.inlet) annotation (Line(
          points={{24,2},{50,2}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
    end SimpleReaction2;

    model SimpleReaction22 "The simple chemical reaction A+B<->C+D with equilibrium [C]*[D]/([A]*[B]) = 2, where [A] is molar concentration of A in water"
      import Chemical;

       extends Modelica.Icons.Example;

      constant Real Kb(unit="kg/mol") = 2
        "Molarity based dissociation constant of the reaction with one more reactant";

      constant Real Kx(unit="1") = Kb*55.508
        "Mole fraction based dissociation constant of the reaction with one more reactant in the pure water";

      constant Modelica.Units.SI.Temperature T_25degC=298.15 "Temperature";
      constant Real R = Modelica.Constants.R "Gas constant";

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Chemical.Boundaries.Substance A(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,2},{-14,22}})));
      Chemical.Processes.Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx)),
            Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-R*T_25degC*log(Kx))},
        nS=2,
        nP=2) annotation (Placement(transformation(extent={{4,-8},{24,12}})));
      Chemical.Boundaries.Substance B(use_mass_start=false, amountOfSubstance_start=0.1) annotation (Placement(transformation(extent={{-34,-24},{-14,-4}})));
      Chemical.Boundaries.Substance C(
        use_mass_start=false,
        amountOfSubstance_start=0.1,
        useInlet=true) annotation (Placement(transformation(extent={{48,-8},{68,12}})));

      Chemical.Boundaries.Substance D(
        use_mass_start=false,
        amountOfSubstance_start=0.1,
        useInlet=true) annotation (Placement(transformation(extent={{44,-34},{64,-14}})));
      inner DropOfCommons dropOfCommons(L=1)    annotation (Placement(transformation(extent={{52,58},{72,78}})));
    equation
      connect(A.solution, solution.solution) annotation (Line(
          points={{-30,2},{-30,-90},{60,-90},{60,-98}},
          color={127,127,0}));
      connect(C.solution, solution.solution) annotation (Line(points={{52,-8},{66,-8},{66,-90},{60,-90},{60,-98}},
                                                 color={127,127,0}));
      connect(B.solution, solution.solution) annotation (Line(points={{-30,-24},
            {-30,-90},{60,-90},{60,-98}},color={127,127,0}));

      connect(B.outlet, reaction2_1.substrates[1]) annotation (Line(
          points={{-14,-14},{-4,-14},{-4,1.75},{4,1.75}},
          color={158,66,200},
          thickness=0.5));
      connect(A.outlet, reaction2_1.substrates[2]) annotation (Line(
          points={{-14,12},{-4,12},{-4,2.25},{4,2.25}},
          color={158,66,200},
          thickness=0.5));
      connect(D.solution, solution.solution) annotation (Line(points={{48,-34},{60,-34},{60,-98}}, color={127,127,0}));
      connect(reaction2_1.products[1], C.inlet) annotation (Line(
          points={{24,1.75},{36,1.75},{36,2},{48,2}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[2], D.inlet) annotation (Line(
          points={{24,2.25},{34,2.25},{34,-24},{44,-24}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Simple reaction demonstrating equilibria between substance A, B, and substance C, mixed in one solution. Observe the molar concentration (A.c) and molar fraction. Note, that molar fractions (A.x and B.x and C.x) are always summed to 1 for the whole solution.</p>
</html>"),
        experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
    end SimpleReaction22;

    model HeatingOfWater "Heating of 1 kg water"
      extends Modelica.Icons.Example;

      Chemical.Solution solution(useMechanicPorts=true, useThermalPort=true) annotation (Placement(transformation(extent={{-98,-100},{102,100}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{-86,-72},{-66,-52}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
      Chemical.Boundaries.Substance liquidWater(
        substanceData=Chemical.SubstancesOld.Water_liquid(),
        use_mass_start=true,
        mass_start=1) annotation (Placement(transformation(extent={{22,-28},{42,-8}})));
      inner Modelica.Fluid.System system(T_ambient=298.15)
        annotation (Placement(transformation(extent={{60,50},{80,70}})));
    equation
      connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
          points={{-66,-62},{-58,-62},{-58,-102}},
          color={191,0,0}));
    connect(fixed1.flange, solution.bottom) annotation (Line(
        points={{-18,-84},{2,-84},{2,-102}},
        color={0,127,0}));
      connect(solution.solution, liquidWater.solution) annotation (Line(points={{62,-98},{26,-98},{26,-28}},
                                          color={127,127,0}));
      annotation (experiment(StopTime=1),
      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Heating of solution by one degree, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T) for temperature change.</p>
</html>"));
    end HeatingOfWater;

    model HeatingOfAlcohol "Heating of 50% ethanol"
      extends Modelica.Icons.Example;

      Chemical.Solution solution(useMechanicPorts=true, useThermalPort=true) annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=4180)
        annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));
      Chemical.Boundaries.Substance Ethanol(
        redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
        substanceData=Chemical.SubstancesOld.Ethanol_liquid(),
        use_mass_start=true,
        mass_start=(55.508/2)*0.04607) annotation (Placement(transformation(extent={{18,-8},{38,12}})));

      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,-94},{-8,-74}})));
      Chemical.Boundaries.Substance liquidWater(
        substanceData=Chemical.SubstancesOld.Water_liquid(),
        use_mass_start=true,
        mass_start=1/2) annotation (Placement(transformation(extent={{-50,-8},{-30,12}})));
      inner Modelica.Fluid.System system(T_ambient=298.15)
        annotation (Placement(transformation(extent={{62,42},{82,62}})));
    equation
      connect(fixedHeatFlow.port, solution.heatPort) annotation (Line(
          points={{-66,-66},{-60,-66},{-60,-102}},
          color={191,0,0}));
    connect(solution.solution, Ethanol.solution) annotation (Line(
        points={{60,-98},{60,-34},{22,-34},{22,-8}},
        color={127,127,0}));
    connect(solution.bottom, fixed1.flange) annotation (Line(
        points={{0,-102},{0,-84},{-18,-84}},
        color={0,127,0}));
      connect(solution.solution, liquidWater.solution) annotation (Line(points={{
              60,-98},{60,-34},{-46,-34},{-46,-8}}, color={127,127,0}));
      annotation (experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Heating of solution of water and ethanol, using standard HeatPort from Modelica Standard Library.</p>
<p>Observe Solution.T (or H2O.Solution.T or Ethanol.Solution.T) for temperature change. Note, that we can heat all substances in solution at once and the results would differ from HeatingOfWater.</p>
</html>"));
    end HeatingOfAlcohol;

    model ExothermicReaction
      "Exothermic reaction in ideally thermal isolated solution and in constant temperature conditions"
      import Chemical;

       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.MolarEnergy ReactionEnthalpy=-55000;

      Chemical.Solution thermal_isolated_solution(useMechanicPorts=true, ConstantTemperature=false)
        annotation (Placement(transformation(extent={{-100,-100},{98,-6}})));
      Chemical.Boundaries.Substance A(use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,-60},{-20,-40}})));
      Chemical.Processes.Reaction reaction2_2(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
      Chemical.Boundaries.Substance B(
        use_mass_start=false,
        amountOfSubstance_start=0.1,
        useInlet=true) annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

      Chemical.Solution solution_at_constant_temperature(useMechanicPorts=true, useThermalPort=true)
        annotation (Placement(transformation(extent={{-100,0},{98,94}})));
      Chemical.Boundaries.Substance A1(use_mass_start=false, amountOfSubstance_start=0.9) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Chemical.Processes.Reaction reaction2_1(
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfH=ReactionEnthalpy)},
        nS=1,
        nP=1) annotation (Placement(transformation(extent={{-8,40},{12,60}})));
      Chemical.Boundaries.Substance B1(
        use_mass_start=false,
        amountOfSubstance_start=0.1,
        useInlet=true) annotation (Placement(transformation(extent={{20,40},{40,60}})));

      //  Modelica.SIunits.HeatFlowRate q
      //    "Heat flow to environment to reach constant temperature";
      Modelica.Units.SI.Temperature t
        "Temperature if the solution is ideally thermal isolated from environment";
      Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{20,4},{40,24}})));
      Chemical.Boundaries.Substance H2O1(
        substanceData=Chemical.SubstancesOld.Water_liquid(),
        use_mass_start=true,
        mass_start=1,
        initAmount=Chemical.Utilities.Types.InitializationMethods.state,
        amountOfSubstance_start=55) annotation (Placement(transformation(extent={{26,-94},{46,-74}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-28,4},{-8,24}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed2
        annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
      inner Modelica.Fluid.System system(T_ambient=298.15)
        annotation (Placement(transformation(extent={{56,64},{76,84}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
            298.15)
        annotation (Placement(transformation(extent={{-88,26},{-68,46}})));
    equation
      //  q = fixedTemperature.port.Q_flow;
      t = thermal_isolated_solution.solution.T;

      connect(B.solution, thermal_isolated_solution.solution) annotation (Line(
          points={{24,-60},{24,-64},{58.4,-64},{58.4,-99.06}},
          color={127,127,0}));
      connect(A.solution, thermal_isolated_solution.solution) annotation (Line(
            points={{-36,-60},{-36,-64},{58.4,-64},{58.4,-99.06}},
                                                             color={127,127,0}));
      connect(B1.solution, solution_at_constant_temperature.solution) annotation (
          Line(
          points={{24,40},{24,34},{58.4,34},{58.4,0.94}},
          color={127,127,0}));
      connect(A1.solution, solution_at_constant_temperature.solution) annotation (
          Line(points={{-36,40},{-36,34},{58.4,34},{58.4,0.94}},
                                                          color={127,127,0}));
    connect(solution_at_constant_temperature.solution, H2O.solution)
      annotation (Line(
        points={{58.4,0.94},{24,0.94},{24,4}},
        color={127,127,0}));
    connect(thermal_isolated_solution.solution, H2O1.solution) annotation (Line(
        points={{58.4,-99.06},{30,-99.06},{30,-94}},
        color={127,127,0}));
    connect(solution_at_constant_temperature.bottom, fixed1.flange) annotation (
       Line(
        points={{-1,-0.94},{0,-0.94},{0,14},{-18,14}},
        color={0,127,0}));
    connect(thermal_isolated_solution.bottom, fixed2.flange) annotation (Line(
        points={{-1,-100.94},{-1,-86},{-16,-86}},
        color={0,127,0}));
      connect(solution_at_constant_temperature.heatPort, fixedTemperature.port)
        annotation (Line(points={{-60.4,-0.94},{-60.4,36},{-68,36}}, color={191,0,
              0}));
      connect(A1.outlet, reaction2_1.substrates[1]) annotation (Line(
          points={{-20,50},{-8,50}},
          color={158,66,200},
          thickness=0.5));
      connect(A.outlet, reaction2_2.substrates[1]) annotation (Line(
          points={{-20,-50},{-8,-50}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_1.products[1], B1.inlet) annotation (Line(
          points={{12,50},{20,50}},
          color={158,66,200},
          thickness=0.5));
      connect(reaction2_2.products[1], B.inlet) annotation (Line(
          points={{12,-50},{20,-50}},
          color={158,66,200},
          thickness=0.5));
      annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>Demonstration of exotermic reaction with perfect cooling (i.e. connected fixed temperature to the HeatPort) and thermally insulated (HetPort unconnected). See solution_(...).T</p>
</html>"),
        experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}})));
    end ExothermicReaction;

    model HydrogenCombustion "Hydrogen combustion in piston"
      import Chemical;
      extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.Volume V=0.001 "Initial volume";
     // parameter Modelica.SIunits.Pressure p=100000 "Initial pressure";
      parameter Modelica.Units.SI.Temperature T=298.15
        "Initial temperature";

      parameter Modelica.Units.SI.Area A=0.01 "Cross area of cylinder";

      //p*V=n*R*T
      // parameter Modelica.SIunits.AmountOfSubstance n=p*V/(Modelica.Constants.R*T)
      //   "Initial amount of substances in sulution";
      Chemical.Solution idealGas(
        SurfaceArea=A,
        useMechanicPorts=true,
        useThermalPort=true,
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas)
                                                               annotation (Placement(transformation(extent={{-108,-50},{-8,50}})));
                       // AmbientPressure=p)
      //  volume_start=V,
      Chemical.Boundaries.Substance H2_gas(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        substanceData=Chemical.SubstancesOld.Hydrogen_gas(),
        use_mass_start=false,
        amountOfSubstance_start=26) annotation (Placement(transformation(extent={{-98,-26},{-78,-6}})));
      Chemical.Boundaries.Substance O2_gas(
        substanceData=Chemical.SubstancesOld.Oxygen_gas(),
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=13) annotation (Placement(transformation(extent={{-100,10},{-80,30}})));
      Chemical.Boundaries.Substance H2O_gas(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=1,
        useInlet=true) annotation (Placement(transformation(extent={{-34,-8},{-14,12}})));
      Chemical.Processes.FastReaction fastReaction(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas "Ideal Gas",
        p={2},
        s={1,2},
        productsSubstanceData={Chemical.SubstancesOld.Water_gas()},
        nS=2,
        nP=1) annotation (Placement(transformation(extent={{-68,-8},{-48,12}})));
      Modelica.Mechanics.Translational.Components.Spring spring(c=1e6) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-58,64})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=200000)
        annotation (Placement(transformation(extent={{-94,-80},{-74,-60}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature coolerTemperature(T=298.15)
        annotation (Placement(transformation(extent={{-18,-80},{-38,-60}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-58,78})));
      Modelica.Mechanics.Translational.Components.Fixed fixed1
        annotation (Placement(transformation(extent={{-68,-66},{-48,-46}})));

      inner DropOfCommons dropOfCommons         annotation (Placement(transformation(extent={{52,62},{72,82}})));
    equation
    connect(H2_gas.solution, idealGas.solution) annotation (Line(
        points={{-94,-26},{-28,-26},{-28,-49}},
        color={127,127,0}));
    connect(O2_gas.solution, idealGas.solution) annotation (Line(
        points={{-96,10},{-102,10},{-102,-26},{-28,-26},{-28,-49}},
        color={127,127,0}));
    connect(H2O_gas.solution, idealGas.solution) annotation (Line(
        points={{-30,-8},{-30,-26},{-28,-26},{-28,-49}},
        color={127,127,0}));
      connect(idealGas.surfaceFlange, spring.flange_a) annotation (Line(
          points={{-58,50},{-58,54}},
          color={0,127,0}));
      connect(idealGas.heatPort, thermalConductor.port_a) annotation (Line(
          points={{-88,-51},{-88,-56},{-106,-56},{-106,-70},{-94,-70}},
          color={191,0,0}));
      connect(thermalConductor.port_b, coolerTemperature.port) annotation (Line(
          points={{-74,-70},{-38,-70}},
          color={191,0,0}));
      connect(fixed.flange, spring.flange_b) annotation (Line(
          points={{-58,78},{-58,74}},
          color={0,127,0}));
    connect(idealGas.bottom, fixed1.flange) annotation (Line(
        points={{-58,-51},{-58,-56}},
        color={0,127,0}));
      connect(O2_gas.outlet, fastReaction.substrates[1])
        annotation (Line(
          points={{-80,20},{-74,20},{-74,1.75},{-68,1.75}},
          color={158,66,200},
          thickness=0.5));
      connect(H2_gas.outlet, fastReaction.substrates[2])
        annotation (Line(
          points={{-78,-16},{-74,-16},{-74,2.25},{-68,2.25}},
          color={158,66,200},
          thickness=0.5));
      connect(fastReaction.products[1], H2O_gas.inlet) annotation (Line(
          points={{-48,2},{-34,2}},
          color={158,66,200},
          thickness=0.5));
      annotation ( experiment(StopTime=0.0009, __Dymola_Algorithm="Dassl"),
                                           Documentation(info="<html>
<p>The gaseous reaction of hydrogen combustion: </p>
<table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
<td><p align=\"center\"><b>2 H<sub>2</sub> + O<sub>2</sub> &lt;-&gt; 2 H<sub>2</sub>O</b></p></td>
<td><p>(1)</p></td>
</tr>
</table>
<p><br>This reaction generates a large amount of energy which can be used for mechanical or thermal purposes. </p>
<p>Building this model using the Chemical library components is easy. First, we drag and drop the library class &lsquo;Components.Solution&rsquo; into the diagram of our new model, labeled &lsquo;idealGas&rsquo; in Figure 4. In parameter dialog of this solution we check &ldquo;useThermalPorts&rdquo; and &ldquo;useMechanicsPorts&rdquo; to enable the thermal and mechanical interface. In the same dialog we need to set the area of the piston (e.g., 1 dm<sup>2</sup>), where the pressure provides the force of the green mechanical port of the uppermost side. The next parameter is the ambient external pressure surrounding the system (e.g., 1 bar). All three chemical substances of the reaction (1) can be added by dragging and dropping the library class &lsquo;Components.Substance&rsquo;. Because this model uses gases, the state of matter must be changed to some gas, such as the ideal gas prepared as &lsquo;Interfaces.IdealGas&rsquo;. The substance data must be selected to define the appropriate substances such as &lsquo;Hydrogen_gas&rsquo;, &lsquo;.Oxygen_gas&rsquo; and &lsquo;.Water_gas&rsquo; in package &lsquo;Examples.Substances&rsquo;. In addition, the initial amounts of substances can be prepared for the ideal solution of hydrogen and oxygen gases at a ratio 2:1 to attain the chemical equation above, with the expectation that at the end of the burning process, only water vapor would be presented. Therefore, the initial values of H<sub>2</sub> particles could be set to 26 mmol and of O<sub>2</sub> particles as 13 mmol. All substances must be connected with the &lsquo;idealGas&rsquo; using the blue colored solution port situated on the bottom side of each substance and solution. Then, the chemical reaction is inserted into the diagram of this model as library class &lsquo;Components.Reaction&rsquo;, and it is set to two substrates (nS=2) with stoichiometry s={2,1} and one product with stoichiometry p={2} to represent the reaction (3). The substances are then connected using violet colored substance connectors with appropriate indexes: H<sub>2</sub> to substrates[1], O<sub>2</sub> to substrates[2] and H<sub>2</sub>O to products[1]. At this point, the model is prepared to simulate the conditions of an unconnected heat port and an unconnected mechanical port. This simulation reaches the theoretical ideal of thermally isolated (zero heat flow from/to the solution) and isobaric (zero force generated on piston) conditions. </p>
<p><br><img src=\"modelica://Chemical/Resources/Images/Examples/HydrogenBurning.png\"/></p>
<p><font style=\"color: #222222; \">Mueller, M. A., Kim, T. J., Yetter, R. A., &amp; Dryer, F. L. (1999). Flow reactor studies and kinetic modeling of the H2/O2 reaction.&nbsp;<i>International Journal of Chemical Kinetics</i>,&nbsp;<i>31</i>(2), 113-125.</font></p>
<p><br>However, in the real world, there is always some thermal energy flow from the solution, and this cooling process can be connected using the thermal connector of the Modelica Standard Library 3.2.1. For example, the simple thermal conductor of thermal conductance 2W/K at a constant temperature environment of 25&deg;C is represented in the model. The mechanical power of the engine can be connected to the robust mechanical model. However, in our example we selected only a very strong mechanical spring with a spring constant of 10<sup>6</sup> N/m to stop the motion of the piston in order to generate the pressure. This standard spring component is situated above the solution in the model diagram. The results of this experiment are shown in Figure 1. </p>
</html>"),
        Diagram(coordinateSystem(extent={{-120,-100},{120,100}})));
    end HydrogenCombustion;

    model WaterVaporization "Evaporation of water"
      import Chemical;
       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.Temperature T_start=273.15
        "Initial temperature";

      Chemical.Solution liquid(temperature_start=T_start, useThermalPort=true) annotation (Placement(transformation(extent={{-98,-98},{-6,-8}})));

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Solution gas(
        temperature_start=T_start,
        useThermalPort=true,
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas)
                                                               annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Chemical.Boundaries.Substance H2O_gaseuous(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        useInlet=true) annotation (Placement(transformation(extent={{8,50},{28,70}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{56,60},{76,80}})));
      Chemical.Boundaries.Substance otherSubstances(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        substanceData=Chemical.SubstancesOld.Oxygen_gas(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{2,28},{22,48}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(
        G=1e6) annotation (Placement(transformation(extent={{44,-12},{64,8}})));
      Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-28,-62},{-48,-42}})));
      inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
        annotation (Placement(transformation(extent={{72,-74},{92,-54}})));
      Chemical.Processes.GasVolatility gasVolatility(
        redeclare package stateOut = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.Water_gas(),
        k_forward=10,
        redeclare package stateOfMatterOut = Chemical.Interfaces.IdealGas "Ideal Gas")
        annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=180,
            origin={-78,26})));
      inner DropOfCommons dropOfCommons(L=1)
        annotation (Placement(transformation(extent={{-92,72},{-72,92}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G=1e6)
               annotation (Placement(transformation(extent={{44,-40},{64,-20}})));
    equation

      connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
          points={{27.6,6.9},{12,6.9},{12,50}},
          color={127,127,0}));
    connect(fixedTemperature.T, clock.y) annotation (Line(
        points={{96,8},{104,8},{104,70},{77,70}},
        color={0,0,127},
        smooth=Smooth.Bezier));
    connect(gas.solution, otherSubstances.solution) annotation (Line(
        points={{27.6,6.9},{6,6.9},{6,28}},
        color={127,127,0}));
    connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
        points={{74,8},{72,8},{72,-2},{64,-2}},
        color={191,0,0}));
    connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
        points={{-27.6,5.1},{-28,5.1},{-28,-2},{44,-2}},
        color={191,0,0}));
      connect(liquid.solution, liquidWater.solution) annotation (Line(points={{-24.4,-97.1},{-24.4,-96.55},{-32,-96.55},{-32,-62}},
                                                             color={127,127,0}));
      connect(liquidWater.outlet,gasVolatility. inlet) annotation (Line(
          points={{-48,-52},{-88,-52},{-88,26}},
          color={158,66,200},
          thickness=0.5));
      connect(thermalConductor1.port_b, fixedTemperature.port) annotation (Line(points={{64,-30},{68,-30},{68,-2},{72,-2},{72,8},{74,8}}, color={191,0,0}));
      connect(thermalConductor1.port_a, liquid.heatPort) annotation (Line(points={{44,-30},{0,-30},{0,-102},{-79.6,-102},{-79.6,-98.9}}, color={191,0,0}));
      connect(gasVolatility.solution, gas.solution) annotation (Line(points={{-67.8,21.8},{-52,21.8},{-52,6},{27.6,6},{27.6,6.9}}, color={127,127,0}));
      connect(gasVolatility.outlet, H2O_gaseuous.inlet)
        annotation (Line(
          points={{-68,26},{-52,26},{-52,60},{8,60}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        experiment(
          StopTime=400,
          Tolerance=1e-08,
          __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>Demonstraiton of water vaporization between two solutions - liquid and gaseous. The temperature is increased in time to illustrate, how the vaporization rate rises in higher temperatures. See liquid.T and liquid.Volume, compared to gas.T and gas.volume.</p>
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterVaporization;

    model WaterSublimation "Sublimation of water"
      import Chemical;
       extends Modelica.Icons.Example;

      parameter Modelica.Units.SI.Temperature T_start=273.15 - 50
        "Initial temperature";

      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Solution gas(
        temperature_start=T_start,
        useThermalPort=true,
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        BasePressure=600) annotation (Placement(transformation(extent={{-46,6},{46,96}})));
                                    /*volume_start(
        displayUnit="l") = 0.001, */
      Chemical.Boundaries.Substance H2O_gaseuous(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=0.001,
        useInlet=true) annotation (Placement(transformation(extent={{-20,56},{0,76}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
                                                           fixedTemperature
                  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          origin={84,8})));
      Modelica.Blocks.Sources.ContinuousClock clock(offset=1*T_start)
        annotation (Placement(transformation(extent={{62,36},{82,56}})));
      Chemical.Boundaries.Substance otherSubstances(
        substanceData=Chemical.SubstancesOld.Oxygen_gas(),
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{-4,36},{16,56}})));
      Chemical.Solution solid(
        temperature_start=T_start,
        BasePressure=600,
        useThermalPort=true) annotation (Placement(transformation(extent={{-50,-100},{42,-10}})));
      Chemical.Boundaries.Substance H2O_solid(
        substanceData=Chemical.SubstancesOld.Water_IceIh(),
        use_mass_start=false,
        amountOfSubstance_start=55.508) "Solid water" annotation (Placement(transformation(extent={{10,-52},{-10,-32}})));
      Chemical.Processes.GasVolatility gasVolatility(
        redeclare package stateOut = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.Water_gas(),
        redeclare package stateOfMatterOut = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-66,26})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=1e6)
        annotation (Placement(transformation(extent={{48,-8},{68,12}})));
      inner DropOfCommons dropOfCommons
        annotation (Placement(transformation(extent={{-88,74},{-68,94}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G=1e6)
        annotation (Placement(transformation(extent={{50,-32},{70,-12}})));
    equation

      connect(gas.solution, H2O_gaseuous.solution) annotation (Line(
          points={{27.6,6.9},{27.6,0},{-50,0},{-50,56},{-16,56}},
          color={127,127,0}));
    connect(fixedTemperature.T, clock.y) annotation (Line(
        points={{96,8},{98,8},{98,46},{83,46}},
        color={0,0,127},
        smooth=Smooth.Bezier));
    connect(gas.solution, otherSubstances.solution) annotation (Line(
        points={{27.6,6.9},{24,6.9},{24,6},{20,6},{20,36},{0,36}},
        color={127,127,0}));
      connect(solid.solution, H2O_solid.solution) annotation (Line(
          points={{23.6,-99.1},{23.6,-104},{-50,-104},{-50,-52},{6,-52}},
          color={127,127,0}));
      connect(fixedTemperature.port, thermalConductor.port_b) annotation (Line(
          points={{74,8},{72,8},{72,2},{68,2}},
          color={191,0,0}));
      connect(gas.heatPort, thermalConductor.port_a) annotation (Line(
          points={{-27.6,5.1},{-28,5.1},{-28,2},{48,2}},
          color={191,0,0}));
      connect(thermalConductor1.port_a, solid.heatPort) annotation (Line(points={{50,-22},{44,-22},{44,-90},{-31.6,-90},{-31.6,-100.9}}, color={191,0,0}));
      connect(thermalConductor1.port_b, fixedTemperature.port) annotation (Line(points={{70,-22},{72,-22},{72,8},{74,8}}, color={191,0,0}));
      connect(H2O_solid.outlet, gasVolatility.inlet) annotation (Line(
          points={{-10,-42},{-76,-42},{-76,26}},
          color={158,66,200},
          thickness=0.5));
      connect(gasVolatility.outlet, H2O_gaseuous.inlet) annotation (Line(
          points={{-56,26},{-56,66},{-20,66}},
          color={200,66,175},
          thickness=0.5));
      connect(gas.solution, gasVolatility.solution) annotation (Line(points={{27.6,6.9},{27.6,12},{-55.8,12},{-55.8,21.8}}, color={127,127,0}));
      annotation (
        experiment(StopTime=49.7, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>Demonstraiton of water sublimation between two solutions - solid and gaseous. The temperature is increased in time to illustrate, how the sublimation rate rises in higher temperatures. See solid.T and solid.Volume, compared to gas.T and gas.volume. Note, that the liquid phase is omitted here.</p>
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end WaterSublimation;

    model GasSolubility_NIST "Dissolution of gases in liquids"
      import Chemical;
       extends Modelica.Icons.Example;

      Chemical.Solution water_solution_25degC(temperature_start=298.15) annotation (Placement(transformation(extent={{-160,-78},{-68,12}})));
                                          //(amountOfSolution_start=52.3)
      Chemical.Solution water_solution_37degC(temperature_start=310.15) annotation (Placement(transformation(extent={{-52,-80},{42,12}})));
                                       //(amountOfSolution_start=39.7)
      Chemical.Processes.GasSolubility CO2_dissolutionP(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{-138,42},{-118,62}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Boundaries.Substance CO2_25(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        useInlet=true) "Free dissolved CO2 in water at 25 degC" annotation (Placement(transformation(extent={{-130,-28},{-150,-8}})));

      Chemical.Processes.GasSolubility O2_dissolutionP(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.Oxygen_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{-100,40},{-80,60}})));

      Chemical.Boundaries.ExternalIdealGas O2_g_25(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure(displayUnit="mmHg") = 12665.626804425)
        annotation (Placement(transformation(extent={{-114,74},{-94,94}})));
      Chemical.Boundaries.Substance O2_25(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) "Free dissolved O2 in water at 25 degC" annotation (Placement(transformation(extent={{-94,-26},{-114,-6}})));

      Chemical.Processes.GasSolubility CO2_dissolutionE(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{-24,40},{-4,60}})));

      Chemical.Boundaries.ExternalIdealGas CO2_g_25(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure(displayUnit="mmHg") =
          5332.8954966) annotation (Placement(transformation(extent={{-154,74},{-134,94}})));

      Chemical.Boundaries.Substance CO2_37(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        useInlet=true) "Free dissolved CO2 in water at 37degC" annotation (Placement(transformation(extent={{-22,-34},{-42,-14}})));

      Chemical.Processes.GasSolubility O2_dissolutionE_NIST(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.Oxygen_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{18,42},{38,62}})));
      Chemical.Boundaries.Substance O2_37(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) "Free dissolved O2 in water at 37degC" annotation (Placement(transformation(extent={{18,-34},{-2,-14}})));

      Chemical.Boundaries.Substance water_25(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-100,-68},{-80,-48}})));
      Chemical.Boundaries.Substance water_37(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{8,-70},{28,-50}})));
      Chemical.Boundaries.ExternalIdealGas CO2_g_37(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure(displayUnit="mmHg") =
          5332.8954966) annotation (Placement(transformation(extent={{-44,68},{-24,88}})));
      Chemical.Boundaries.ExternalIdealGas O2_g_37(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure(displayUnit="mmHg") = 12665.626804425)
        annotation (Placement(transformation(extent={{-6,68},{14,88}})));
      Solution water_solution_37degC1(temperature_start=273.15) annotation (Placement(transformation(extent={{66,-80},{160,12}})));
      Chemical.Processes.GasSolubility CO2_dissolutionE1(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{92,44},{112,64}})));
      Chemical.Boundaries.Substance CO2_0(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        useInlet=true) "Free dissolved CO2 in water at 0degC" annotation (Placement(transformation(extent={{96,-34},{76,-14}})));

      Chemical.Processes.GasSolubility O2_dissolutionE_NIST1(
        redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        substanceDataOut=Chemical.SubstancesOld.Oxygen_aqueous(),
        redeclare package stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas",
        k_forward=1) annotation (Placement(transformation(extent={{134,42},{154,62}})));
      Chemical.Boundaries.Substance O2_0(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) "Free dissolved O2 in water at 0degC" annotation (Placement(transformation(extent={{136,-34},{116,-14}})));

      Chemical.Boundaries.Substance water_0(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{126,-70},{146,-50}})));
      Chemical.Boundaries.ExternalIdealGas CO2_g_0(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure(displayUnit="mmHg") = 5332.8954966)
        annotation (Placement(transformation(extent={{74,68},{94,88}})));
      Chemical.Boundaries.ExternalIdealGas O2_g_0(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure(displayUnit="mmHg") = 12665.626804425)
        annotation (Placement(transformation(extent={{112,68},{132,88}})));
      inner Modelica.Fluid.System system(p_ambient=100000)
        annotation (Placement(transformation(extent={{-70,-98},{-50,-78}})));
      Real kH_CO2_25, kH_O2_25;
      inner DropOfCommons dropOfCommons
        annotation (Placement(transformation(extent={{38,76},{58,96}})));
    equation

      kH_CO2_25 = CO2_25.c / CO2_g_25.substance.x;
      kH_O2_25 = O2_25.c / O2_g_25.substance.x;
    //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;
    //  kH_CO2_25 = CO2_25.x / CO2_g_25.x;

      connect(CO2_25.solution, water_solution_25degC.solution) annotation (Line(
            points={{-134,-28},{-134,-77.1},{-86.4,-77.1}},
                                                          color={127,127,0}));
      connect(O2_25.solution, water_solution_25degC.solution) annotation (Line(
            points={{-98,-26},{-98,-77.1},{-86.4,-77.1}}, color={127,127,0}));
      connect(CO2_37.solution, water_solution_37degC.solution) annotation (Line(
            points={{-26,-34},{-26,-79.08},{23.2,-79.08}},
                                                         color={127,127,0}));
      connect(O2_37.solution, water_solution_37degC.solution) annotation (Line(
            points={{14,-34},{14,-79.08},{23.2,-79.08}}, color={127,127,0}));
      connect(water_25.solution, water_solution_25degC.solution) annotation (Line(
            points={{-96,-68},{-96,-77.1},{-86.4,-77.1}}, color={127,127,0}));
      connect(water_37.solution, water_solution_37degC.solution) annotation (Line(
            points={{12,-70},{12,-79.08},{23.2,-79.08}}, color={127,127,0}));
      connect(CO2_0.solution, water_solution_37degC1.solution) annotation (Line(
            points={{92,-34},{92,-79.08},{141.2,-79.08}},   color={127,127,0}));
      connect(O2_0.solution, water_solution_37degC1.solution) annotation (Line(
            points={{132,-34},{132,-79.08},{141.2,-79.08}}, color={127,127,0}));
      connect(water_0.solution, water_solution_37degC1.solution) annotation (Line(
            points={{130,-70},{130,-79.08},{141.2,-79.08}}, color={127,127,0}));
      connect(CO2_g_25.outlet, CO2_dissolutionP.inlet) annotation (Line(
          points={{-134,84},{-138,84},{-138,52}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_g_25.outlet, O2_dissolutionP.inlet) annotation (Line(
          points={{-94,84},{-100,84},{-100,50}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2_g_37.outlet, CO2_dissolutionE.inlet) annotation (Line(
          points={{-24,78},{-24,50}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_g_37.outlet, O2_dissolutionE_NIST.inlet) annotation (Line(
          points={{14,78},{18,78},{18,52}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2_g_0.outlet, CO2_dissolutionE1.inlet) annotation (Line(
          points={{94,78},{92,78},{92,54}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_g_0.outlet, O2_dissolutionE_NIST1.inlet) annotation (Line(
          points={{132,78},{132,52},{134,52}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2_g_25.solution, water_solution_25degC.solution)
        annotation (Line(points={{-150,74},{-150,70},{-156,70},{-156,-78},{-122,-78},{-122,-77.1},{-86.4,-77.1}}, color={127,127,0}));
      connect(O2_g_25.solution, water_solution_25degC.solution)
        annotation (Line(points={{-110,74},{-110,68},{-156,68},{-156,-78},{-122,-78},{-122,-77.1},{-86.4,-77.1}}, color={127,127,0}));
      connect(CO2_g_37.solution, water_solution_37degC.solution) annotation (Line(points={{-40,68},{-46,68},{-46,-79.08},{23.2,-79.08}}, color={127,127,0}));
      connect(O2_g_37.solution, water_solution_37degC.solution) annotation (Line(points={{-2,68},{-46,68},{-46,-79.08},{23.2,-79.08}}, color={127,127,0}));
      connect(CO2_g_0.solution, water_solution_37degC1.solution) annotation (Line(points={{78,68},{70,68},{70,-79.08},{141.2,-79.08}}, color={127,127,0}));
      connect(O2_g_0.solution, water_solution_37degC1.solution) annotation (Line(points={{116,68},{70,68},{70,-79.08},{141.2,-79.08}}, color={127,127,0}));
      connect(CO2_dissolutionP.outlet, CO2_25.inlet) annotation (Line(
          points={{-118,52},{-118,-18},{-130,-18}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_dissolutionP.outlet, O2_25.inlet) annotation (Line(
          points={{-80,50},{-80,-16},{-94,-16}},
          color={200,66,175},
          thickness=0.5));
      connect(CO2_dissolutionE.outlet, CO2_37.inlet) annotation (Line(
          points={{-4,50},{-4,-24},{-22,-24}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_dissolutionE_NIST.outlet, O2_37.inlet) annotation (Line(
          points={{38,52},{38,-24},{18,-24}},
          color={200,66,175},
          thickness=0.5));
      connect(CO2_dissolutionE1.outlet, CO2_0.inlet)
        annotation (Line(
          points={{112,54},{112,-26},{96,-26},{96,-24}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_dissolutionE_NIST1.outlet, O2_0.inlet)
        annotation (Line(
          points={{154,52},{154,-26},{136,-26},{136,-24}},
          color={200,66,175},
          thickness=0.5));
      connect(CO2_dissolutionP.solution, water_solution_25degC.solution)
        annotation (Line(points={{-117.8,47.8},{-117.8,34},{-76,34},{-76,-78},{-86.4,-78},{-86.4,-77.1}}, color={127,127,0}));
      connect(water_solution_25degC.solution, O2_dissolutionP.solution)
        annotation (Line(points={{-86.4,-77.1},{-86.4,-78},{-76,-78},{-76,45.8},{-79.8,45.8}}, color={127,127,0}));
      connect(water_solution_37degC.solution, CO2_dissolutionE.solution)
        annotation (Line(points={{23.2,-79.08},{23.2,-80},{40,-80},{40,36},{-2,36},{-2,45.8},{-3.8,45.8}}, color={127,127,0}));
      connect(water_solution_37degC.solution, O2_dissolutionE_NIST.solution)
        annotation (Line(points={{23.2,-79.08},{23.2,-80},{40,-80},{40,47.8},{38.2,47.8}}, color={127,127,0}));
      connect(CO2_dissolutionE1.solution, water_solution_37degC1.solution)
        annotation (Line(points={{112.2,49.8},{112.2,34},{168,34},{168,-79.08},{141.2,-79.08}}, color={127,127,0}));
      connect(O2_dissolutionE_NIST1.solution, water_solution_37degC1.solution)
        annotation (Line(points={{154.2,47.8},{168,47.8},{168,-79.08},{141.2,-79.08}}, color={127,127,0}));
      annotation (
        experiment(StopTime=1, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>Demonstration of CO2 and O2 dissolution in pure water as described by NIST Henry&apos;s law data at 25degC.</p>
<p>Recalculation from Henry&apos;s law constants ( https://webbook.nist.gov/cgi/inchi?ID=C124389&amp;Mask=10#Solubility ):</p>
<p>CO2:</p>
<p>kH = 0.035 mol/(kg.bar) at 25degC</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 1.87 mmol/kg at 25degC</li>
</ul>
<p><br>Other temperatures:</p>
<p>NIST constant ... 2400 K ( https://webbook.nist.gov/cgi/inchi?ID=C124389&amp;Mask=10#Solubility )</p>
<p>0degC ... 0.035*exp(2400*(1/273.15 - 1/298.15) = 0.073 mol/(kg.bar)</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 3.9 mmol/kg at 0degC</li>
</ul>
<p>37degC ... 0.035*exp(2400*(1/310.15 - 1/298.15) = 0.026 mol/(kg.bar)</p>
<ul>
<li>pCO2 = 40 mmHg =&gt; dissolved CO2 .. 1.4 mmol/kg at 37degC</li>
</ul>
<p><br>O2:</p>
<p>kH = 0.0013 mol/(kg.bar) at 25degC</p>
<ul>
<li>pO2 = 95 mmHg (0.126656 bar) =&gt; dissolved O2 .. 0.165 mmol/kg at 25degC</li>
</ul>
<p><br>Other temperatures:</p>
<p>NIST constant ... 1500 K ( https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&amp;Mask=10 )</p>
<p>0degC ... 0.0013*exp(1500*(1/273.15 - 1/298.15) = 0.0021 mol/(kg.bar)</p>
<ul>
<li>pO2 = 95 mmHg =&gt; dissolved O2 .. 0.26 mmol/kg at 0degC</li>
</ul>
<p>37degC ... 0.0013*exp(1500*(1/273.15 - 1/298.15) = 0.0011 mol/(kg.bar)</p>
<ul>
<li>pO2 = 95 mmHg =&gt; dissolved O2 .. 0.14 mmol/kg at 37degC</li>
</ul>
</html>", revisions="<html>
<p><i>2015-2019</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),
        Diagram(coordinateSystem(extent={{-160,-100},{160,100}})));
    end GasSolubility_NIST;

    model GasSolubilityInBlood "Dissolution of gases in liquids"
      import Chemical;

       extends Modelica.Icons.Example;

      Chemical.Solution blood_plasma annotation (Placement(transformation(extent={{-100,-76},{-8,14}})));
                                          //(amountOfSolution_start=52.3)
      Chemical.Solution red_cells annotation (Placement(transformation(extent={{8,-78},{102,14}})));
                                       //(amountOfSolution_start=39.7)
      Chemical.Processes.GasSolubility CO2_dissolutionP(redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.CarbonDioxide_aqueous()) annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
      //  kH_T0(displayUnit="(mol/kg H2O)/bar at 25degC,101325Pa")= 0.00062064026806947,
      Chemical.Boundaries.Substance CO2_unbound_plasma(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        useInlet=true) "Free dissolved CO2 in blood plasma" annotation (Placement(transformation(extent={{-70,-26},{-90,-6}})));

      Chemical.Processes.GasSolubility O2_dissolutionP(redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.Oxygen_aqueous()) annotation (Placement(transformation(extent={{-34,44},{-14,64}})));

      Chemical.Boundaries.ExternalIdealGas O2_g_n1(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure=12665.626804425)
        annotation (Placement(transformation(extent={{22,78},{42,98}})));
      Chemical.Boundaries.Substance O2_unbound_plasma(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) "Free dissolved O2 in blood plasma" annotation (Placement(transformation(extent={{-30,-28},{-50,-8}})));

      Chemical.Processes.GasSolubility CO2_dissolutionE(redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.CarbonDioxide_aqueous()) annotation (Placement(transformation(extent={{36,44},{56,64}})));

      Chemical.Boundaries.ExternalIdealGas CO2_g_n2(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure(displayUnit="mmHg") =
          5332.8954966) annotation (Placement(transformation(extent={{-58,78},{-38,98}})));

      Chemical.Boundaries.Substance CO2_unbound_erythrocyte(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.001,
        useInlet=true) "Free dissolved CO2 in red cells" annotation (Placement(transformation(extent={{38,-34},{18,-14}})));

      Chemical.Processes.GasSolubility O2_dissolutionE_NIST(redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.Oxygen_aqueous()) annotation (Placement(transformation(extent={{78,44},{98,64}})));
      Chemical.Boundaries.Substance O2_unbound_erythrocyte_NIST(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) "Free dissolved O2 in red cells" annotation (Placement(transformation(extent={{78,-32},{58,-12}})));

      Chemical.Boundaries.Substance water_plasma(
        substanceData=Chemical.SubstancesOld.Water_liquid(),
        use_mass_start=true,
        mass_start=0.82) annotation (Placement(transformation(extent={{-40,-66},{-20,-46}})));
      Chemical.Boundaries.Substance water_erythrocyte(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=0.66)
        annotation (Placement(transformation(extent={{72,-68},{92,-48}})));
      inner Modelica.Fluid.System system(p_ambient(displayUnit="mmHg")=
          101325.0144354, T_ambient=310.15)
        annotation (Placement(transformation(extent={{-10,-96},{10,-76}})));
      Chemical.Boundaries.Substance other_plasma(
        use_mass_start=true,
        mass_start=0.18,
        substanceData=Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=1/0.627))
        annotation (Placement(transformation(extent={{-70,-60},{-50,-40}})));
      Chemical.Boundaries.Substance other_erythrocyte(mass_start=0.34, substanceData=Chemical.Interfaces.Incompressible.SubstanceDataParameters(MolarWeight=1.7,
            density=1200)) annotation (Placement(transformation(extent={{38,-68},{58,-48}})));
    equation

    connect(CO2_unbound_plasma.solution, blood_plasma.solution) annotation (
        Line(
        points={{-74,-26},{-74,-75.1},{-26.4,-75.1}},
        color={127,127,0}));
    connect(O2_unbound_plasma.solution, blood_plasma.solution) annotation (Line(
        points={{-34,-28},{-34,-75.1},{-26.4,-75.1}},
        color={127,127,0}));
    connect(CO2_unbound_erythrocyte.solution, red_cells.solution) annotation (
        Line(
        points={{34,-34},{34,-77.08},{83.2,-77.08}},
        color={127,127,0}));
    connect(O2_unbound_erythrocyte_NIST.solution, red_cells.solution)
      annotation (Line(
        points={{74,-32},{74,-77.08},{83.2,-77.08}},
        color={127,127,0}));
      connect(water_plasma.solution, blood_plasma.solution) annotation (Line(points=
             {{-36,-66},{-36,-75.1},{-26.4,-75.1}}, color={127,127,0}));
      connect(water_erythrocyte.solution, red_cells.solution) annotation (Line(points={{76,-68},{76,-77.08},{83.2,-77.08}}, color={127,127,0}));
      connect(other_plasma.solution, blood_plasma.solution) annotation (Line(points={{-66,-60},{-66,-75.1},{-26.4,-75.1}}, color={127,127,0}));
      connect(CO2_g_n2.solution, blood_plasma.solution) annotation (Line(points={{-54,78},{-96,78},{-96,-75.1},{-26.4,-75.1}}, color={127,127,0}));
      connect(O2_g_n1.solution, red_cells.solution) annotation (Line(points={{26,78},{12,78},{12,-77.08},{83.2,-77.08}}, color={127,127,0}));
      connect(CO2_g_n2.outlet, CO2_dissolutionP.inlet)
        annotation (Line(
          points={{-38,88},{-28,88},{-28,70},{-78,70},{-78,54}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2_g_n2.outlet, CO2_dissolutionE.inlet)
        annotation (Line(
          points={{-38,88},{-28,88},{-28,70},{36,70},{36,54}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_g_n1.outlet, O2_dissolutionP.inlet)
        annotation (Line(
          points={{42,88},{52,88},{52,74},{-34,74},{-34,54}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_g_n1.outlet, O2_dissolutionE_NIST.inlet)
        annotation (Line(
          points={{42,88},{52,88},{52,74},{78,74},{78,54}},
          color={158,66,200},
          thickness=0.5));
      connect(CO2_dissolutionP.outlet, CO2_unbound_plasma.inlet)
        annotation (Line(
          points={{-58,54},{-58,-16},{-70,-16}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_dissolutionP.outlet, O2_unbound_plasma.inlet)
        annotation (Line(
          points={{-14,54},{-14,-18},{-30,-18}},
          color={200,66,175},
          thickness=0.5));
      connect(CO2_dissolutionE.outlet, CO2_unbound_erythrocyte.inlet)
        annotation (Line(
          points={{56,54},{56,-24},{38,-24}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_dissolutionE_NIST.outlet, O2_unbound_erythrocyte_NIST.inlet)
        annotation (Line(
          points={{98,54},{98,-22},{78,-22},{78,-22}},
          color={200,66,175},
          thickness=0.5));
      connect(blood_plasma.solution, CO2_dissolutionP.solution)
        annotation (Line(points={{-26.4,-75.1},{-96,-75.1},{-96,38},{-46,38},{-46,49.8},{-57.8,49.8}}, color={127,127,0}));
      connect(O2_dissolutionP.solution, blood_plasma.solution)
        annotation (Line(points={{-13.8,49.8},{-6,49.8},{-6,38},{-96,38},{-96,-75.1},{-26.4,-75.1}}, color={127,127,0}));
      connect(CO2_dissolutionE.solution, red_cells.solution)
        annotation (Line(points={{56.2,49.8},{64,49.8},{64,38},{12,38},{12,-77.08},{83.2,-77.08}}, color={127,127,0}));
      connect(O2_dissolutionE_NIST.solution, red_cells.solution)
        annotation (Line(points={{98.2,49.8},{106,49.8},{106,38},{12,38},{12,-77.08},{83.2,-77.08}}, color={127,127,0}));
      connect(O2_unbound_erythrocyte_NIST.solution, water_erythrocyte.solution) annotation (Line(points={{74,-32},{74,-68},{76,-68}}, color={127,127,0}));
      connect(other_erythrocyte.solution, red_cells.solution) annotation (Line(points={{42,-68},{42,-78},{40,-78},{40,-77.08},{83.2,-77.08}}, color={127,127,0}));
      annotation (
        experiment(StopTime=1, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>Demonstration of different blood gases solubility in erythrocytes and in plasma. The difference is governed by various amount of other substances in the solution. </p>
<p>Please note, that the total content of CO2 and O2 in blood plasma and erythrocytes must be determined by including bicarbonate and hemoglobin connected amounts. </p>
</html>", revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end GasSolubilityInBlood;

    model EnzymeKinetics "Basic enzyme kinetics"
      import Chemical;

      extends Modelica.Icons.Example;

      Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

      //The huge negative Gibbs energy of the product will make the second reaction almost irreversible (e.g. K=exp(50))
      Chemical.Boundaries.Substance P(
        useInlet=true,
        useSolution=true,
        mass_start=1e-8,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{72,-12},{92,8}})));

      Chemical.Boundaries.Substance S(use_mass_start=false, amountOfSubstance_start=100) annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));

      parameter Modelica.Units.SI.AmountOfSubstance tE=1
        "Total amount of enzyme";
         parameter Real k_cat(
        unit="mol/s",
        displayUnit="mol/min")=1
        "Forward rate of second reaction";
      constant Modelica.Units.SI.Concentration Km=0.1
        "Michaelis constant = substrate concentration at rate of half Vmax";

      parameter Modelica.Units.SI.MolarFlowRate Vmax=1e-5*k_cat
        "Maximal molar flow";

      Chemical.Boundaries.Substance ES(
        useInlet=true,
        mass_start=tE/2,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{-8,-10},{12,10}})));
      Chemical.Boundaries.Substance E(
        useInlet=true,
        mass_start=tE/2,
        amountOfSubstance_start=tE/2) annotation (Placement(transformation(extent={{10,36},{-10,56}})));
      Chemical.Processes.Reaction chemicalReaction(
        k_forward=1,
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*log(2/Km))},
        nP=1,
        nS=2) annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      Chemical.Processes.ForwardReaction chemicalReaction1(
        k_forward=k_cat,
        productsSubstanceData={Chemical.Interfaces.Incompressible.SubstanceDataParameters(DfG=-Modelica.Constants.R*298.15*50),
            Chemical.Interfaces.Incompressible.SubstanceDataParameters()},
        nS=1,
        nP=2) annotation (Placement(transformation(extent={{24,-12},{44,8}})));

      Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{42,-80},{62,-60}})));
      inner DropOfCommons dropOfCommons      annotation (Placement(transformation(extent={{68,70},{88,90}})));
    equation
      //Michaelis-Menton: v=((E.q_out.conc + ES.q_out.conc)*k_cat)*S.concentration/(Km+S.concentration);
      connect(E.solution, solution.solution) annotation (Line(
          points={{6,36},{-8,36},{-8,-98},{60,-98}},
          color={127,127,0}));
      connect(ES.solution, solution.solution)
        annotation (Line(points={{-4,-10},{-4,-98},{60,-98}},         color={127,127,0}));

      connect(S.solution, solution.solution) annotation (Line(
          points={{-88,-14},{-88,-56},{-8,-56},{-8,-98},{60,-98}},
          color={127,127,0}));
      connect(P.solution, solution.solution) annotation (Line(
          points={{76,-12},{76,-98},{60,-98}},
          color={127,127,0}));
      connect(liquidWater.solution, solution.solution) annotation (Line(points={{
              46,-80},{46,-98},{60,-98}}, color={127,127,0}));
      connect(ES.outlet, chemicalReaction1.substrates[1]) annotation (Line(
          points={{12,0},{16,0},{16,-2},{24,-2}},
          color={158,66,200},
          thickness=0.5));
      connect(chemicalReaction1.products[1], P.inlet)
        annotation (Line(
          points={{44,-2.25},{58,-2.25},{58,-2},{72,-2}},
          color={200,66,175},
          thickness=0.5));
      connect(chemicalReaction1.products[2], E.inlet)
        annotation (Line(
          points={{44,-1.75},{60,-1.75},{60,46},{10,46}},
          color={200,66,175},
          thickness=0.5));
      connect(chemicalReaction.products[1], ES.inlet) annotation (Line(
          points={{-20,0},{-8,0}},
          color={158,66,200},
          thickness=0.5));
      connect(S.outlet, chemicalReaction.substrates[1])
        annotation (Line(
          points={{-72,-4},{-66,-4},{-66,-0.25},{-40,-0.25}},
          color={158,66,200},
          thickness=0.5));
      connect(E.outlet, chemicalReaction.substrates[2])
        annotation (Line(
          points={{-10,46},{-52,46},{-52,0},{-40,0},{-40,0.25}},
          color={158,66,200},
          thickness=0.5));
          annotation ( Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",     info="<html>
<p>Be carefull, the assumption for Michaelis-Menton are very strong: </p>
<p>The substrate must be in sufficiently high concentration and the product must be in very low concentration to reach almost all enzyme in enzyme-substrate complex all time. ([S] &gt;&gt; Km) &amp;&amp; ([P] &lt;&lt; K2)</p>
<p><br>To recalculate the enzyme kinetics from Michaelis-Menton parameters Km, tE a k_cat is selected the same half-rate of the reaction defined as:</p>
<p>E = ES = tE/2 .. the amount of free enzyme is the same as the amount of enzyme-substrate complexes</p>
<p>S = Km .. the amount of substrate is Km</p>
<p>r = Vmax/2 = tE*k_cat / 2 .. the rate of reaction is the half of maximal rate</p>
<p><br>Conversions of molar concentration to mole fraction (MM is molar mass of the solvent in solution -&gt; 55.508 kg/mol for water):</p>
<p>x(Km) = Km/MM</p>
<p>x(tE) = tE/MM</p>
<p>xS = S/MM = Km/MM</p>
<p><br>The new kinetics of the system defined as:</p>
<p>uS&deg; = DfG(S) = 0</p>
<p>uE&deg; = DfG(E) = 0</p>
<p>uES&deg; = <b>DfG(ES) = DfG(S) + DfG(E) - R*T*ln(2/x(Km))</b></p>
<p>from dissociation coeficient of the frist reaction 2/x(Km) = xSE/(xS*xE) = exp((uE&deg; + uS&deg; - uES&deg;)/(RT))</p>
<p>uP&deg; = DfG(P) </p>
<p><br>r = Vmax/2</p>
<p>r = -kC1 * (uES&deg; - uE&deg; - uS&deg; + R*T*ln(xES/(xE*xS) ) = -kC1 * (-R*T*ln(2/x(Km)) + R*T*ln(xS) ) = kC1 * R * T * ln(2)</p>
<p>because xES=xE this time</p>
<p>r = -kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) = -kC2 * (DfG(P) - uES&deg; + R*T*ln(xP) ) = kC2 * (-DfG(P) - R * T * ln(2))</p>
<h4>kC1 = (Vmax/2) / (R * T * ln(2))</h4>
<h4>kC2 = (Vmax/2) / ( -DfG(P) - R * T * ln(2) ) </h4>
<p><br>For example in case of C=AmountOfSolution/(Tau*ActivationPotential) we can rewrite C to ActivationPotential (Be carefull: this energy is not the same as in <a href=\"http://en.wikipedia.org/wiki/Arrhenius_equation\">Arrhenius equation</a> or in Transition State Theory):</p>
<p>ActivationPotential1 = AmountOfSolution/(Tau*(Vmax/2)) * R * T * ln(2) </p>
<p>ActivationPotential2 = AmountOfSolution/(Tau*(Vmax/2)) * ( -DfG(P) - R * T * ln(2) ) </p>
<p><br>where</p>
<p>AmountOfSolution = MM = 55.508 (for water)</p>
<p>Tau = 1 s (just to be physical unit correct)</p>
<p>DfG(P) = -R*T*50 is Gibbs energy of formation of product (setting negative enough makes second reaction almost irreversible)</p>
<h4>The maximum of the new enzyme kinetics</h4>
<p>The enzymatic rate must have a maximum near of Vmax. </p>
<p>The new maximum is a litle higher: Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) ), for example if -uP&deg;/RT = 50, the new maximum is around 1.014*Vmax, where Vmax is the maximum of Michaelis Menten.</p>
<p>The proof:</p>
<p>We want to sutisfied the following inequality:</p>
<p>-kC2 * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) ?=&lt;? Vmax * (1 + 1/( -uP&deg;/(R*T*ln(2)) - 1) )</p>
<p><br>(Vmax/2) * (uP&deg; + uE&deg; - uES&deg; + R*T*ln(xP*xE/xES) ) / ( - uP&deg; - R * T * ln(2) ) ?=&lt;? Vmax*(1 + R*T*ln(2) / ( -uP&deg; - R*T*ln(2)) )</p>
<p>(uP&deg; +<b> </b>R*T*ln(2/x(Km)) + R*T*ln(xP*xE/xES) ) ?=&lt;? 2*( - uP&deg; - R * T * ln(2) ) + 2*R*T*ln(2)</p>
<p>R*T*ln(xP*xE/xES) ?=&lt;? - uP&deg; - R*T*ln(2/x(Km)) </p>
<p>xP*xE/xES ?=&lt;? exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>The equality is the equation of the equilibrium: xP*xE/xES = exp((- uP&deg; - uE&deg; + uES&deg; )/(R*T)) = exp((- uP&deg; - R*T*ln(2/x(Km))/(R*T))</p>
<p>If the equilibrium of the reaction is reached only by forward rate then xP*xE/xES must be less than the dissociation constant.</p>
</html>"),
        experiment(StopTime=5000, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end EnzymeKinetics;

    model WaterElectrolysis "Water electrolysis"
      import Chemical;

      extends Modelica.Icons.Example;
      Chemical.Boundaries.Substance O2_gas(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=0.001,
        useInlet=true) annotation (Placement(transformation(extent={{-12,-6},{8,14}})));

      Chemical.Boundaries.Substance H2_gas(
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
        use_mass_start=false,
        amountOfSubstance_start=0.001,
        useInlet=true) annotation (Placement(transformation(extent={{38,-4},{18,16}})));
      Chemical.Processes.Reaction reaction(
        s={2,4},
        p={2,1,4},
        productsSubstanceData={Chemical.SubstancesOld.Hydrogen_aqueous(),Chemical.SubstancesOld.Oxygen_aqueous(),Chemical.SubstancesOld.Electrone_solid()},
        nS=2,
        nP=3) annotation (Placement(transformation(
            extent={{11,11},{-11,-11}},
            rotation=180,
            origin={-31,-29})));
      Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,-78},{92,30}})));
      Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-80},{-56,28}})));
      Chemical.Solution water(temperature_start=310.15) annotation (Placement(transformation(extent={{-28,-84},{52,-42}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-42,70},{-22,90}})));
      Chemical.Boundaries.ElectronSink electrone annotation (Placement(transformation(extent={{84,-38},{64,-18}})));
      Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{-84,-34},{-64,-14}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
      annotation (Placement(transformation(extent={{-36,38},{-16,58}})));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
      annotation (Placement(transformation(extent={{-66,38},{-46,58}})));
      Chemical.Solution air(redeclare package stateOfMatter =
            Chemical.Interfaces.IdealGas)                                                   annotation (Placement(transformation(extent={{-40,-16},{50,26}})));
      Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
        annotation (Placement(transformation(extent={{18,38},{-2,58}})));
      Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{0,-74},{-20,-54}})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{36,26},{56,46}})));
      Chemical.Boundaries.Substance O2_aq(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001,
        useInlet=true) annotation (Placement(transformation(extent={{4,-68},{24,-48}})));
      Chemical.Boundaries.Substance H2_aq(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 1e-05,
        useInlet=true) annotation (Placement(transformation(extent={{28,-72},{48,-52}})));
      Chemical.Processes.GasVolatility H2_volatility(redeclare package stateOut = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.Hydrogen_gas()) annotation (Placement(transformation(extent={{64,-12},{44,8}})));
      Chemical.Processes.GasVolatility O2_volatility(redeclare package stateOut = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
            Chemical.SubstancesOld.Oxygen_gas()) annotation (Placement(transformation(extent={{-38,-18},{-18,2}})));
    equation
    connect(electrone1.pin,voltageSensor. p) annotation (Line(
        points={{-74,-14.2},{-92,-14.2},{-92,48},{-74,48},{-74,80},{-42,80}},
        color={0,0,255}));
    connect(electrone.pin,voltageSensor. n) annotation (Line(
        points={{74,-18.2},{74,48},{60,48},{60,80},{-22,80}},
        color={0,0,255}));
    connect(electrone.solution,anode. solution) annotation (Line(
        points={{80,-38},{80,-76.92},{85.2,-76.92}},
        color={127,127,0}));
    connect(electrone1.pin,currentSensor. p) annotation (Line(
        points={{-74,-14.2},{-92,-14.2},{-92,48},{-66,48}},
        color={0,0,255}));
    connect(currentSensor.n,resistor. p) annotation (Line(
        points={{-46,48},{-36,48}},
        color={0,0,255}));
    connect(electrone1.solution,cathode. solution) annotation (Line(
        points={{-80,-34},{-80,-66},{-74,-66},{-74,-78.92},{-62.8,-78.92}},
        color={127,127,0}));
      connect(O2_gas.solution, air.solution) annotation (Line(points={{-8,-6},{-8,-20},{32,-20},{32,-15.58}},
                                             color={127,127,0}));
      connect(H2_gas.solution, air.solution) annotation (Line(points={{34,-4},{34,-15.58},{32,-15.58}},
                            color={127,127,0}));
      connect(constantVoltage.p, voltageSensor.n) annotation (Line(points={{18,48},{
              60,48},{60,80},{-22,80}}, color={0,0,255}));
      connect(resistor.n, constantVoltage.n)
        annotation (Line(points={{-16,48},{-2,48}}, color={0,0,255}));
      connect(liquidWater.solution, water.solution) annotation (Line(points={{-4,-74},{-4,-88},{36,-88},{36,-83.58}},
                                               color={127,127,0}));
      connect(ground.p, constantVoltage.p) annotation (Line(points={{46,46},{46,48},{18,48}}, color={0,0,255}));
      connect(liquidWater.outlet, reaction.substrates[1])
        annotation (Line(
          points={{-20,-64},{-46,-64},{-46,-29.275},{-42,-29.275}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone1.outlet, reaction.substrates[2])
        annotation (Line(
          points={{-64,-24},{-48,-24},{-48,-28.725},{-42,-28.725}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_aq.solution, water.solution) annotation (Line(points={{8,-68},{8,-88},{36,-88},{36,-83.58}}, color={127,127,0}));
      connect(H2_aq.solution, water.solution) annotation (Line(points={{32,-72},{32,-78},{36,-78},{36,-83.58}}, color={127,127,0}));
      connect(reaction.products[1], H2_aq.inlet)
        annotation (Line(
          points={{-20,-29.3667},{-20,-30},{28,-30},{28,-62}},
          color={200,66,175},
          thickness=0.5));
      connect(reaction.products[2], O2_aq.inlet)
        annotation (Line(
          points={{-20,-29},{-18,-29},{-18,-30},{4,-30},{4,-58}},
          color={200,66,175},
          thickness=0.5));
      connect(reaction.products[3], electrone.inlet)
        annotation (Line(
          points={{-20,-28.6333},{-18,-28},{64,-28}},
          color={200,66,175},
          thickness=0.5));
      connect(H2_aq.outlet, H2_volatility.inlet) annotation (Line(
          points={{48,-62},{52,-62},{52,-2},{64,-2}},
          color={158,66,200},
          thickness=0.5));
      connect(H2_volatility.outlet, H2_gas.inlet) annotation (Line(
          points={{44,-2},{42,-2},{42,6},{38,6}},
          color={200,66,175},
          thickness=0.5));
      connect(O2_aq.outlet, O2_volatility.inlet) annotation (Line(
          points={{24,-58},{20,-58},{20,-40},{-42,-40},{-42,-8},{-38,-8}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_volatility.solution, air.solution) annotation (Line(points={{-17.8,-12.2},{-18,-12.2},{-18,-16},{32,-16},{32,-15.58}}, color={127,127,0}));
      connect(H2_volatility.solution, air.solution) annotation (Line(points={{43.8,-6.2},{43.8,-20},{32,-20},{32,-15.58}}, color={127,127,0}));
      connect(O2_volatility.outlet, O2_gas.inlet) annotation (Line(
          points={{-18,-8},{-16,-8},{-16,4},{-12,4}},
          color={158,66,200},
          thickness=0.5));
      annotation ( experiment(StopTime=1), Documentation(info="<html>
<p>The water ecectrolysis: </p>
<p><b>2 H<sub>2</sub>O +&nbsp;&nbsp;4 e<sup>-</sup><sub>(catode)</sub>&nbsp;&lt;-&gt;  2 H<sub>2</sub> + O<sub>2</sub>&nbsp;+&nbsp;&nbsp;4 e<sup>-</sup><sub>(anode)</sub>&nbsp;</b></p>
</html>"));
    end WaterElectrolysis;

    model ElectrochemicalCell
      "The electrochemical cell: Pt(s) | H2(g) | H+(aq), Cl-(aq) | AgCl(s) | Ag(s)"
      import Chemical;

     extends Modelica.Icons.Example;

      Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
      Chemical.Boundaries.Substance Ag(
        use_mass_start=false,
        amountOfSubstance_start=1,
        useInlet=true) annotation (Placement(transformation(extent={{-52,-30},{-72,-10}})));

      Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

      Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-60},{38,6}})));

      Chemical.Boundaries.Substance Cl(
        use_mass_start=false,
        useInlet=true,
        initAmount=Chemical.Utilities.Types.InitializationMethods.state,
        amountOfSubstance_start=12.39) annotation (Placement(transformation(extent={{-14,-26},{6,-6}})));

      Chemical.Boundaries.Substance AgCl(
        substanceData=Chemical.SubstancesOld.SilverChloride_solid(),
        use_mass_start=false,
        amountOfSubstance_start=1e-8) annotation (Placement(transformation(extent={{-76,4},{-56,24}})));

      Chemical.Boundaries.Substance H(
        substanceData=Chemical.SubstancesOld.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=12.39) annotation (Placement(transformation(extent={{10,-26},{30,-6}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Processes.FastReaction electrodeReaction(
        s={2,2},
        productsSubstanceData={Chemical.SubstancesOld.Hydrogen_aqueous()},
        nS=2,
        nP=1) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=270,
            origin={52,6})));
      Chemical.Processes.FastReaction electrodeReaction1(
        productsSubstanceData={Chemical.SubstancesOld.Silver_solid(),Chemical.SubstancesOld.Chloride_aqueous()},
        nP=2,
        nS=2) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={-40,0})));

      Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{84,-84},{104,-64}})));
      Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
        annotation (Placement(transformation(extent={{-6,-54},{14,-34}})));
      inner DropOfCommons dropOfCommons      annotation (Placement(transformation(extent={{-80,-82},{-60,-62}})));
      Chemical.Processes.Membrane membrane annotation (Placement(transformation(extent={{-40,-44},{-20,-24}})));
      Chemical.Boundaries.ExternalIdealGas            H2(
        useInlet=true,
        useSolution=false,
        redeclare package gasModel = Chemical.Interfaces.IdealGas "Ideal Gas",
        PartialPressure=100000,
        redeclare package stateOfMatter = Chemical.Interfaces.IdealGas "Ideal Gas")
                  annotation (Placement(transformation(extent={{60,52},{80,72}})));
    equation
      connect(Cl.solution, solution1.solution) annotation (Line(
          points={{-10,-26},{-10,-30},{24.4,-30},{24.4,-59.34}},
          color={127,127,0}));
      connect(H.solution, solution1.solution) annotation (Line(points={{14,-26},{14,-30},{24.4,-30},{24.4,-59.34}},
                                         color={127,127,0}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,-26},{84,-49},{89.2,-49}},
        color={127,127,0}));
    connect(AgCl.solution, cathode.solution) annotation (Line(
        points={{-72,4},{-74,4},{-74,-34},{-68,-34},{-68,-42.84},{-54.4,-42.84}},
        color={127,127,0}));
    connect(Ag.solution, cathode.solution) annotation (Line(
        points={{-56,-30},{-56,-42.84},{-54.4,-42.84}},
        color={158,66,200}));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-6,74},{-96,74},{-96,51.8},{-68,51.8}},
          color={0,0,255}));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{14,74},{92,74},{92,-6.2},{78,-6.2}},
          color={0,0,255}));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{78,-6.2},{92,-6.2},{92,-64},{94,-64}},
          color={0,0,255}));
      connect(liquidWater.solution, solution1.solution) annotation (Line(points={
              {-2,-54},{-2,-59.34},{24.4,-59.34}}, color={127,127,0}));
      connect(electrone1.outlet, electrodeReaction.substrates[1])
        annotation (Line(
          points={{68,-16},{52.25,-16},{52.25,-4}},
          color={158,66,200},
          thickness=0.5));
      connect(H.outlet, electrodeReaction.substrates[2])
        annotation (Line(
          points={{30,-16},{51.75,-16},{51.75,-4}},
          color={158,66,200},
          thickness=0.5));
      connect(Ag.inlet, electrodeReaction1.products[1])
        annotation (Line(
          points={{-52,-20},{-40.25,-20},{-40.25,-10}},
          color={200,66,175},
          thickness=0.5));
      connect(AgCl.outlet, electrodeReaction1.substrates[1])
        annotation (Line(
          points={{-56,14},{-40.25,14},{-40.25,10}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone.outlet, electrodeReaction1.substrates[2])
        annotation (Line(
          points={{-58,42},{-39.75,42},{-39.75,10}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone.solution, cathode.solution) annotation (Line(points={{-74,32},{-86,32},{-86,30},{-96,30},{-96,-42.84},{-54.4,-42.84}}, color={127,127,0}));
      connect(electrodeReaction1.products[2], membrane.inlet)
        annotation (Line(
          points={{-39.75,-10},{-40,-12},{-40,-34}},
          color={158,66,200},
          thickness=0.5));
      connect(membrane.outlet, Cl.inlet)
        annotation (Line(
          points={{-20,-34},{-18,-34},{-18,-16},{-14,-16}},
          color={158,66,200},
          thickness=0.5));
      connect(solution1.solution, membrane.solution)
        annotation (Line(points={{24.4,-59.34},{24.4,-60},{-2,-60},{-2,-64},{-19.8,-64},{-19.8,-38.2}}, color={127,127,0}));
      connect(electrodeReaction.products[1], H2.inlet) annotation (Line(
          points={{52,16},{52,62},{60,62}},
          color={158,66,200},
          thickness=0.5));
      annotation (
      experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end ElectrochemicalCell;

    model LeadAcidBattery
      "The electrochemical cell: PbSO4(s) | Pb(s) | HSO4-(aq) , H+(aq) | PbO2(s) | PbSO4(s) + 2 H2O"
      import Chemical;

     extends Modelica.Icons.Example;

      Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{26,-74},{60,34}})));

      Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-80,-78},{-46,30}})));

      Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-26,-80},{2,20}})));

      Chemical.Boundaries.Substance Pb(
        substanceData=Chemical.SubstancesOld.Lead_solid(),
        use_mass_start=false,
        amountOfSubstance_start=50) annotation (Placement(transformation(extent={{52,-66},{32,-46}})));

      Chemical.Boundaries.Substance HSO4(
        substanceData=Chemical.SubstancesOld.HydrogenSulfate_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1) annotation (Placement(transformation(extent={{4,-70},{-16,-50}})));
      Chemical.Boundaries.Substance PbSO4_(
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mol") = 0.001,
        useInlet=true) annotation (Placement(transformation(extent={{52,-30},{32,-10}})));
      Chemical.Boundaries.Substance H(
        use_mass_start=false,
        amountOfSubstance_start=1,
        useInlet=true) annotation (Placement(transformation(extent={{0,-42},{-20,-22}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-32,72},{-12,92}})));
      Chemical.Processes.FastReactionWithSolutions electrodeReaction(
        s={1,1,3,2},
        p={1,2},
        productsSubstanceData={Chemical.SubstancesOld.LeadSulfate_solid(),Chemical.SubstancesOld.Water_liquid()},
        nP=2,
        nS=4) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-34,-14})));
      Chemical.Processes.FastReactionWithSolutions electrodeReaction1(
        p={1,1,2},
        productsSubstanceData={Chemical.SubstancesOld.LeadSulfate_solid(),Chemical.SubstancesOld.Proton_aqueous(),Chemical.SubstancesOld.Electrone_solid()},
        nP=3,
        nS=2) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={18,-14})));

      Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{-78,-38},{-58,-18}})));
      Chemical.Boundaries.Substance PbO2(
        substanceData=Chemical.SubstancesOld.LeadDioxide_solid(),
        use_mass_start=false,
        amountOfSubstance_start=50) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-60,-58})));
      Chemical.Boundaries.Substance H2O(
        use_mass_start=false,
        amountOfSubstance_start=0.114/0.018015,
        useInlet=true) annotation (Placement(transformation(extent={{-22,-6},{-2,14}})));
      Chemical.Boundaries.Substance PbSO4(
        use_mass_start=false,
        useInlet=true,
        initAmount=Chemical.Utilities.Types.InitializationMethods.steadyState,
        amountOfSubstance_start=0.001) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={-60,6})));

    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{16,30},{36,50}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=1)
      annotation (Placement(transformation(extent={{-14,40},{6,60}})));
    Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
      annotation (Placement(transformation(extent={{-56,40},{-36,60}})));

      //Real density, molality, totalmolality, voltage;
      inner Modelica.Fluid.System system(T_ambient=299.15)
        annotation (Placement(transformation(extent={{62,64},{82,84}})));
      inner DropOfCommons dropOfCommons(L=1e-5)
                                             annotation (Placement(transformation(extent={{70,-80},{90,-60}})));
      Chemical.Topology.SplitterT1 splitterT1 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-26,-48})));
      Chemical.Boundaries.ElectronSink electron annotation (Placement(transformation(extent={{54,10},{34,30}})));
    equation
      /*density = solution1.solution.m/solution1.solution.V;
  totalmolality = solution1.solution.n/((H2O.x*solution1.solution.n)*H2O.substanceDataVar.MolarWeight);
  molality = HSO4.x*solution1.solution.n/((H2O.x*solution1.solution.n)*H2O.substanceDataVar.MolarWeight);
  voltage = voltageSensor.v;
*/
      connect(HSO4.solution, solution1.solution) annotation (Line(
          points={{0,-70},{-6,-70},{-6,-78},{-3.6,-78},{-3.6,-79}},
          color={127,127,0}));
      connect(H.solution, solution1.solution) annotation (Line(points={{-4,-42},{-4,-78},{-3.6,-78},{-3.6,-79}},
                                         color={127,127,0}));
      connect(H2O.solution, solution1.solution) annotation (Line(
          points={{-18,-6},{-18,-79},{-3.6,-79}},
          color={127,127,0}));
    connect(Pb.solution, anode.solution) annotation (Line(
        points={{48,-66},{48,-72.92},{53.2,-72.92}},
        color={127,127,0}));
    connect(PbSO4_.solution, anode.solution) annotation (Line(
        points={{48,-30},{48,-72.92},{53.2,-72.92}},
        color={127,127,0}));
    connect(PbO2.solution, cathode.solution) annotation (Line(
        points={{-66,-68},{-66,-70},{-60,-70},{-60,-76.92},{-52.8,-76.92}},
        color={127,127,0}));
    connect(electrone1.pin, voltageSensor.p) annotation (Line(
        points={{-68,-18.2},{-82,-18.2},{-82,50},{-64,50},{-64,82},{-32,82}},
        color={0,0,255}));
    connect(electrone1.pin, currentSensor.p) annotation (Line(
        points={{-68,-18.2},{-82,-18.2},{-82,50},{-56,50}},
        color={0,0,255}));
    connect(currentSensor.n, resistor.p) annotation (Line(
        points={{-36,50},{-14,50}},
        color={0,0,255}));
    connect(PbSO4.solution, cathode.solution) annotation (Line(
        points={{-54,-4},{-54,-70},{-60,-70},{-60,-76.92},{-52.8,-76.92}},
        color={127,127,0}));

      connect(electrodeReaction.products[1], PbSO4.inlet)
        annotation (Line(
          points={{-34.25,-4},{-36,-4},{-36,6},{-50,6}},
          color={200,66,175},
          thickness=0.5));
      connect(electrodeReaction.products[2], H2O.inlet)
        annotation (Line(
          points={{-33.75,-4},{-33.75,4},{-22,4}},
          color={200,66,175},
          thickness=0.5));
      connect(electrone1.solution, cathode.solution) annotation (Line(points={{-74,-38},{-92,-38},{-92,-72},{-52.8,-72},{-52.8,-76.92}}, color={127,127,0}));
      connect(electron.solution, anode.solution) annotation (Line(points={{50,10},{70,10},{70,-72},{53.2,-72},{53.2,-72.92}}, color={127,127,0}));
      connect(electron.pin, resistor.n) annotation (Line(points={{44,29.8},{44,56},{12,56},{12,50},{6,50}}, color={0,0,255}));
      connect(electron.pin, ground.p) annotation (Line(points={{44,29.8},{44,56},{26,56},{26,50}}, color={0,0,255}));
      connect(electron.pin, voltageSensor.n) annotation (Line(points={{44,29.8},{44,56},{12,56},{12,82},{-12,82}}, color={0,0,255}));
      connect(electrodeReaction1.products[1], PbSO4_.inlet)
        annotation (Line(
          points={{18.3333,-4},{68,-4},{68,-20},{52,-20}},
          color={200,66,175},
          thickness=0.5));
      connect(electrodeReaction1.products[2], H.inlet) annotation (Line(
          points={{18,-4},{8,-4},{8,-32},{0,-32}},
          color={200,66,175},
          thickness=0.5));
      connect(electrodeReaction1.products[3], electron.inlet)
        annotation (Line(
          points={{17.6667,-4},{17.6667,20},{34,20}},
          color={200,66,175},
          thickness=0.5));
      connect(HSO4.outlet, splitterT1.inlet) annotation (Line(
          points={{-16,-60},{-26,-60},{-26,-58}},
          color={158,66,200},
          thickness=0.5));
      connect(Pb.outlet, electrodeReaction1.substrates[1])
        annotation (Line(
          points={{32,-56},{18.25,-56},{18.25,-24}},
          color={158,66,200},
          thickness=0.5));
      connect(PbO2.outlet, electrodeReaction.substrates[1])
        annotation (Line(
          points={{-50,-58},{-42,-58},{-42,-28},{-34.375,-28},{-34.375,-24}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterT1.outletA, electrodeReaction.substrates[2])
        annotation (Line(
          points={{-36,-48},{-36,-24},{-34.125,-24}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterT1.outletB, electrodeReaction1.substrates[2])
        annotation (Line(
          points={{-16,-48},{4,-48},{4,-44},{17.75,-44},{17.75,-24}},
          color={158,66,200},
          thickness=0.5));
      connect(cathode.solution, electrodeReaction.productSolution[1])
        annotation (Line(points={{-52.8,-76.92},{-52.8,-82},{-44,-82},{-44,-4},{-38.25,-4}}, color={127,127,0}));
      connect(solution1.solution, electrodeReaction.productSolution[2])
        annotation (Line(points={{-3.6,-79},{-3.6,-84},{-44,-84},{-44,-4},{-37.75,-4}}, color={127,127,0}));
      connect(anode.solution, electrodeReaction1.productSolution[1])
        annotation (Line(points={{53.2,-72.92},{53.2,-38},{54,-38},{54,-2},{22,-2},{22,-4},{22.3333,-4}}, color={127,127,0}));
      connect(solution1.solution, electrodeReaction1.productSolution[2])
        annotation (Line(points={{-3.6,-79},{-3.6,-84},{20,-84},{20,0},{22,0},{22,-4}}, color={127,127,0}));
      connect(anode.solution, electrodeReaction1.productSolution[3])
        annotation (Line(points={{53.2,-72.92},{53.2,-36},{54,-36},{54,0},{22,0},{22,-4},{21.6667,-4}}, color={127,127,0}));
      connect(H.outlet, electrodeReaction.substrates[3])
        annotation (Line(
          points={{-20,-32},{-36,-32},{-36,-24},{-33.875,-24}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone1.outlet, electrodeReaction.substrates[4])
        annotation (Line(
          points={{-58,-28},{-36,-28},{-36,-24},{-33.625,-24}},
          color={158,66,200},
          thickness=0.5));
      annotation (
      experiment(StopTime=47455, __Dymola_Algorithm="Lsodar"),
                                  Documentation(revisions=
                        "<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",   info="<html>
<p>The lead-acid electrochemical cells are characterized by two chemical reactions:</p>
<table width=100%>
<tr><th>PbO2 + HSO4- + 3 H+ +2 e- &harr; PbSO4 + 2 H2O</th><td>(1)</td></tr>
<tr><th>Pb + HSO4- &harr; PbSO4 + H+ + 2 e-</th><td>(2)</td></tr>
</table>
<p>The building of one cell of a lead-acid battery starts with the definition of three solutions: two for the lead elec-trodes and one for the liquid-acid solution (Figure 1A). This can be done by dragging and dropping the library class &lsquo;Components.Solution&rsquo; into the diagram. We called the first instance &ldquo;cathode&rdquo;, the second &ldquo;solution&rdquo; and the last &ldquo;anode&rdquo;. We set the parameter &lsquo;Electri-calGround&rsquo; as &ldquo;false&rdquo; for all of these solutions in order to attain the possibility of non-zero voltages. Now we can specify the chemical substances inside the chemical solutions. We drag and drop the library class &lsquo;Compo-nents.Substance&rsquo; into the &ldquo;solution&rdquo; as chemical sub-stances (Figure 1B). H2O(liquid), H+(aqueous) and HSO4-(aqueous) representing the liquid aqueous solu-tion of sulfuric acid. PbSO4(solid) and PbO2(solid) are placed in the &ldquo;cathode&rdquo;, representing the elements of the positive electrode. The substances Pb(solid) and aP-bSO4(solid) are placed into the &ldquo;anode&rdquo;, representing the elements of the negative electrode. All of these sub-stances must be given unique names (e.g., &ldquo;PbSO4&rdquo; for the cathode and &ldquo;aPbSO4&rdquo; for the anode), because the Modelica language does not support two instances with the same name in a single class.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/LeadAcidBatterry1.png\"/></p>
<p>Figure 1) The building of one electro-chemical cell of a lead-acid battery in four steps: A) adding chemical solutions, B) adding chemical substances, C) adding electron transfers and D) adding chemical reactions.</p>
<p>As mentioned above, the appropriate substance data for all these substances must be selected as predefined parametric records, e.g., &lsquo;Exam-ples.Substances.Water_liquid&rsquo;, &lsquo;.Lead_solid&rsquo;, &lsquo;.Lead_dioxide_solid&rsquo;, &lsquo;.Lead_sulfate_solid&rsquo;, and so on. The last, very special substance to be included is an electron. This class is called &lsquo;Compo-nents.ElectronTransfer&rsquo; and it must be added in order for each electrode to transfer electron from the chemical reaction to the electric circuit (Figure 1C). Each of these substances must be connected to the appropriate solu-tion using a solution port situated in the bottom of the component&rsquo;s icons to indicate that they are all mixed in the solution. By having all these substances, it is possi-ble to implement the chemical reactions. Dragging and dropping the library class &lsquo;Components.Reaction&rsquo; for both chemical reactions, and setting their parameters as an appropriate number of reactants, products and stoi-chiometry, allows the connection of each substance with the reaction, as expressed in reaction (1) and reaction (2). This setting can be done using the parameter dialog of the cathode chemical reaction (1) as there are four types of substrates (nS=4) with stoichiometric coeffi-cients: one for the first and second reactant, three for the third reactant and two for the fourth reactant (s={1,1,3,2}). There are also two types of products (nP=2) with stoichiometry: one for PbSO4 and two for water (p={1,2}), following the chemical scheme of the first chemical reaction above. After setting the number of reactants and products, it is possible to connect the substances with reactions. Each instance of reaction has an array of connectors for substrates and an array of con-nectors for products; the user must be very careful to connect each element of these arrays in the same order as defined by stoichiometric coefficients. This means that, for example, the water must be connected in index 2 to products of the first chemical reaction, because we had already selected the order of products by setting the array of stoichiometric coefficients in reaction (1). The chemical reaction (2) must be set analogically as nS=2, nP=3, p={1,1,2} with connections of substance ports of Pb to substrate[1], HSO4- to substrate[2], PbSO4 to prod-uct[1], H+ to product[2] and e- to product[3], as repre-sented in Figure 1D.</p>
<p>The electrochemical cell has already been imple-mented at this stage. However, the simulation requires the initial state of substances, which for the fully charged battery means that almost all elements of the cathode are PbO2 and almost all elements of the anode are Pb. In this state, the sulfuric acid can be concen-trated, which increases the effectiveness of the electro-chemical cell. To set this state, it is possible to just dou-ble-click on PbO2 and Pb and set the amount, e.g., 1mol. To set the pure concentrated sulfuric acid we can also set the amount of SO4- and H+ as 1mol. This fully charged ideal state is ready to simulate when it is con-nected to the electric ground via one of the electric ports of the one electron transfer component.</p>
<p>These batteries can be connected to any electrical cir-cuit that is slowly discharging. For example, if we only connect the simple electric resistance of 1 Ohm as ex-pressed in Figure 1D, then the simulation of the dis-charging process over 13 hours and 45 minutes gives the results of electric current and electric potential, as can be seen in Figure 2. The exchange of the resistor with a voltage source can simulate the charging process for a discharged cell.</p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/LeadAcidBatterry2.png\"/></p>
<p>Figure 2) Discharging simulation of the lead-acid battery cell from Figure 2D, with the initial amount of substances as described in the text.</p>
</html>"),
        __Dymola_experimentSetupOutput);
    end LeadAcidBattery;

    package AcidBase

      model WaterSelfIonization "H2O  <->  OH-   +   H+ "
        import Chemical;
          extends Modelica.Icons.Example;

        Chemical.Solution solution annotation (Placement(transformation(extent={{-72,2},{76,96}})));
        Chemical.Solution solution1 annotation (Placement(transformation(extent={{-76,-98},{72,-4}})));
        Chemical.Boundaries.Substance H3O(
          use_mass_start=false,
          amountOfSubstance_start=1e-7,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={30,70})));

        Chemical.Boundaries.Substance OH(
          use_mass_start=false,
          amountOfSubstance_start=1e-7,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={30,26})));

        Chemical.Boundaries.Substance H2O(mass_start=1, substanceData=Chemical.SubstancesOld.Water_liquid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-30,46})));
        Chemical.Processes.Reaction waterDissociation(
          s={2},
          productsSubstanceData={Chemical.SubstancesOld.Hydronium_aqueous(),Chemical.SubstancesOld.Hydroxide_aqueous()},
          nS=1,
          nP=2) annotation (Placement(transformation(extent={{-12,36},{8,56}})));
              Real pH, pH3O;
        Chemical.Boundaries.Substance H_(
          use_mass_start=false,
          amountOfSubstance_start=1e-7,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-30})));

        Chemical.Boundaries.Substance OH_(
          use_mass_start=false,
          amountOfSubstance_start=1e-7,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={28,-76})));

        Chemical.Boundaries.Substance H2O_(mass_start=1, substanceData=Chemical.SubstancesOld.Water_liquid())
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-32,-56})));
        Chemical.Processes.Reaction waterDissociation_(
          productsSubstanceData={Chemical.SubstancesOld.Proton_aqueous(),Chemical.SubstancesOld.Hydroxide_aqueous()},
          nS=1,
          nP=2) annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));

        inner Modelica.Fluid.System system(p_ambient=100000, T_ambient=298.15)
          annotation (Placement(transformation(extent={{-90,50},{-70,70}})));
      equation
        pH3O = -log10( H3O.substance.a);

        pH = -log10( H_.substance.a);

        connect(H2O.solution, solution.solution) annotation (Line(
            points={{-36,36},{46.4,36},{46.4,2.94}},
            color={127,127,0}));
        connect(OH.solution, solution.solution) annotation (Line(
            points={{24,16},{24,2.94},{46.4,2.94}},
            color={127,127,0}));
        connect(H3O.solution, solution.solution) annotation (Line(
            points={{24,60},{24,2.94},{46.4,2.94}},
            color={127,127,0}));
        connect(H2O_.solution, solution1.solution) annotation (Line(
            points={{-38,-66},{42.4,-66},{42.4,-97.06}},
            color={127,127,0}));
        connect(OH_.solution, solution1.solution) annotation (Line(
            points={{22,-86},{22,-97.06},{42.4,-97.06}},
            color={127,127,0}));
        connect(H_.solution, solution1.solution) annotation (Line(
            points={{22,-40},{22,-97.06},{42.4,-97.06}},
            color={127,127,0}));
        connect(H2O.outlet, waterDissociation.substrates[1]) annotation (Line(
            points={{-20,46},{-12,46}},
            color={158,66,200},
            thickness=0.5));
        connect(waterDissociation.products[1], H3O.inlet) annotation (Line(
            points={{8,45.75},{8,70},{20,70}},
            color={200,66,175},
            thickness=0.5));
        connect(waterDissociation.products[2], OH.inlet) annotation (Line(
            points={{8,46.25},{8,26},{20,26}},
            color={200,66,175},
            thickness=0.5));
        connect(H2O_.outlet, waterDissociation_.substrates[1]) annotation (Line(
            points={{-22,-56},{-14,-56}},
            color={158,66,200},
            thickness=0.5));
        connect(waterDissociation_.products[1], H_.inlet) annotation (Line(
            points={{6,-56.25},{6,-30},{18,-30}},
            color={200,66,175},
            thickness=0.5));
        connect(waterDissociation_.products[2], OH_.inlet)
          annotation (Line(
            points={{6,-55.75},{6,-76},{18,-76}},
            color={200,66,175},
            thickness=0.5));
        annotation ( Documentation(info="<html>
<p>Self-ionization of water.</p>
<p>Ions difference (SID) in water causes the acidity/basicity, where pH = -log10(aH+). An activity of hydrogen ions aH+ is approximated with concentration (mol/l) of the oxonium cations H3O+.</p>
<pre><b>plotExpression(apply(-log10(WaterSelfIonization.H3O.solute)),&nbsp;false,&nbsp;&quot;pH&quot;,&nbsp;1);</b></pre>
<p><br>The titration slope der(pH)/der(SID)=1.48e+6 1/(mol/L) at pH=7.4.</p>
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1));
      end WaterSelfIonization;

      model CarbonDioxideInWater "CO2 as alone acid-base buffer"
        import Chemical;
          extends Modelica.Icons.Example;
          parameter Real KC=1e4;

        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,46}})));
        Chemical.Boundaries.Substance HCO3(
          use_mass_start=false,
          amountOfSubstance_start(displayUnit="mmol") = 1e-08,
          useInlet=true) annotation (Placement(transformation(extent={{-16,-4},{4,16}})));

        Chemical.Processes.FastReaction HendersonHasselbalch(
          TC=dropOfCommons.L,
          productsSubstanceData={Chemical.SubstancesOld.Bicarbonate_aqueous(),Chemical.SubstancesOld.Proton_aqueous()},
          kC=1e-4,
          nS=2,
          nP=2) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-46,-6},{-26,14}})));
        Chemical.Boundaries.ExternalIdealGas CO2_gas(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure(displayUnit="mmHg") =
            5332.8954966) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-76,82})));
        Chemical.Boundaries.Substance H(
          use_mass_start=false,
          amountOfSubstance_start=2e-7,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={10,-30})));

                                              /*(C=2400, kH_T0(
        displayUnit="(mmol/l)/kPa at 25degC") = 0.81805576878885)*/
        Real pH;

        Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-82,-60},{-62,-40}})));
        inner Modelica.Fluid.System system(T_ambient=310.15)
          annotation (Placement(transformation(extent={{48,64},{68,84}})));

        inner Chemical.DropOfCommons dropOfCommons(L=1e-3) annotation (Placement(transformation(extent={{66,12},{86,32}})));
        Chemical.Processes.GasSolubility gasSolubility(redeclare package stateIn = Chemical.Interfaces.IdealGas "Ideal Gas", substanceDataOut=
              Chemical.SubstancesOld.CarbonDioxide_aqueous()) annotation (Placement(transformation(extent={{-70,42},{-50,62}})));
      equation
        pH = -log10( H.substance.a);

        connect(HCO3.solution, solution.solution) annotation (Line(points={{-12,-4},
              {-12,-98.54},{60,-98.54}},color={127,127,0}));
        connect(liquidWater.solution, solution.solution) annotation (Line(points={{-78,-60},{-104,-60},{-104,-104},{60,-104},{60,-98.54}},
                                                     color={127,127,0}));
        connect(H.solution, solution.solution) annotation (Line(points={{4,
                -40},{-12,-40},{-12,-98.54},{60,-98.54}}, color={127,127,0}));
        connect(CO2_gas.solution, solution.solution) annotation (Line(points={{-86,88},{-106,88},{-106,-60},{-104,-60},{-104,-104},{60,-104},{60,-98.54}},
                                                                                                                                 color={127,127,0}));
        connect(HendersonHasselbalch.products[1], HCO3.inlet)
          annotation (Line(
            points={{-26,3.75},{-26,6},{-16,6}},
            color={200,66,175},
            thickness=0.5));
        connect(HendersonHasselbalch.products[2], H.inlet)
          annotation (Line(
            points={{-26,4.25},{-26,0},{-18,0},{-18,-30},{0,-30}},
            color={200,66,175},
            thickness=0.5));
        connect(liquidWater.outlet, HendersonHasselbalch.substrates[1])
          annotation (Line(
            points={{-62,-50},{-54,-50},{-54,3.75},{-46,3.75}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2_gas.outlet, gasSolubility.inlet)
          annotation (Line(
            points={{-76,72},{-76,60},{-78,60},{-78,52},{-70,52}},
            color={158,66,200},
            thickness=0.5));
        connect(gasSolubility.outlet, HendersonHasselbalch.substrates[2])
          annotation (Line(
            points={{-50,52},{-46,52},{-46,4.25}},
            color={158,66,200},
            thickness=0.5));
        connect(gasSolubility.solution, solution.solution)
          annotation (Line(points={{-49.8,47.8},{-49.8,40},{-48,40},{-48,-104},{60,-104},{60,-98.54}}, color={127,127,0}));
        annotation ( Documentation(info="<html>
<p>CO2 solution in water without any other acid-base buffers.</p>
<pre><b>plotExpression(apply(-log10(CarbonDioxideInWater.H3O.solute)),&nbsp;false,&nbsp;&quot;pH&quot;,&nbsp;1);</b></pre>
<p><br>Please note, that OH- (and CO3^-2) can be neglected from electroneutrality calculation, because of very small concentrations (in physiological pH) anyway. </p>
<p>And if SID&gt;0 then also H3O+ can be also neglected from electroneutrality, because only bicarbonate anions HCO3- (or CO3^-2) are needed there to balance the electroneutrality.</p>
<p><br>The partial pressure of CO2 in gas are input parameter. Outputs are an amount of free dissolved CO2 in liquid and an amount of HCO3-.</p>
<p><br>The titration slope der(pH)/der(SID)=17.5 1/(mol/L) at pH=7.4 and pCO2=40 mmHg.</p>
<p><br>Molar heat of formation (aqueous):</p>
<p>CO2:        -413.5 kJ/mol  (gas: -393.5 kJ/mol )</p>
<p>H2O:        -285.8 kJ/mol</p>
<p>HCO3-:        -692.0 kJ/mol</p>
<p>CO3^-2:        -677.1 kJ/mol</p>
<p><br>Enthalphy of reaction H2O + CO2 &lt;-&gt; HCO3- + H+  :         7.3 kJ/mol</p>
<p>Enthalphy of reaction HCO3- &lt;-&gt; CO3^-2 + H+  :        14.9 kJ/mol</p>
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.02));
      end CarbonDioxideInWater;

      model Phosphate
        import Chemical;
          extends Modelica.Icons.Example;
          parameter Real KC=1e-5;
        Chemical.Solution solution annotation (Placement(transformation(extent={{-98,-100},{100,100}})));

        Chemical.Boundaries.Substance H(
          use_mass_start=false,
          useInlet=true,
          initAmount=Chemical.Utilities.Types.InitializationMethods.state,
          amountOfSubstance_start=10^(-7.4)) "hydrogen ions activity" annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={16,-10})));

        Chemical.Boundaries.Substance H3PO4(
          substanceData=Chemical.SubstancesOld.PhosphoricAcid_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start(displayUnit="mol") = 1e-4) annotation (Placement(transformation(extent={{-92,-58},{-72,-38}})));
        Chemical.Boundaries.Substance H2PO4(
          use_mass_start=false,
          useInlet=true,
          initAmount=Chemical.Utilities.Types.InitializationMethods.state,
          amountOfSubstance_start=0.025) annotation (Placement(transformation(extent={{-38,-56},{-18,-36}})));
        Chemical.Boundaries.Substance HPO4(
          use_mass_start=false,
          useInlet=true,
          initAmount=Chemical.Utilities.Types.InitializationMethods.state,
          amountOfSubstance_start=0.006) annotation (Placement(transformation(extent={{16,-56},{36,-36}})));
        Chemical.Boundaries.Substance PO4(
          use_mass_start=false,
          useInlet=true,
          initAmount=Chemical.Utilities.Types.InitializationMethods.state,
          amountOfSubstance_start=1e-08) annotation (Placement(transformation(extent={{78,-72},{98,-52}})));

        Chemical.Processes.FastReaction chemicalReaction(
          productsSubstanceData={Chemical.SubstancesOld.DihydrogenPhosphate_aqueous(),Chemical.SubstancesOld.Proton_aqueous()},
          kC=KC,
          nP=2,
          nS=1) "10^(-1.915 + 3)" annotation (Placement(transformation(extent={{-64,-58},{-44,-38}})));
        Chemical.Processes.FastReaction chemicalReaction1(
          productsSubstanceData={Chemical.SubstancesOld.HydrogenPhosphate_aqueous(),Chemical.SubstancesOld.Proton_aqueous()},
          kC=KC,
          nS=1,
          nP=2) "10^(-6.66 + 3)" annotation (Placement(transformation(extent={{-14,-58},{6,-38}})));
        Chemical.Processes.FastReaction chemicalReaction2(
          productsSubstanceData={Chemical.SubstancesOld.Phosphate_aqueous(),Chemical.SubstancesOld.Proton_aqueous()},
          kC=1e-8,
          nS=1,
          nP=2) "10^(-11.78 + 3)" annotation (Placement(transformation(extent={{44,-58},{64,-38}})));

        Chemical.Boundaries.Substance H2O(use_mass_start=true, amountOfSubstance_start=1/0.018015)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={70,28})));
        Real pH "acidity";
        Chemical.Topology.JunctionN          junctionN(N=3) annotation (Placement(transformation(extent={{-26,-18},{-6,2}})));
      equation
        pH = -log10( H.substance.a);
        connect(H3PO4.solution, solution.solution) annotation (Line(
            points={{-88,-58},{-88,-106},{60.4,-106},{60.4,-98}}));
        connect(H2PO4.solution, solution.solution) annotation (Line(points={{-34,-56},{-34,-88},{60.4,-88},{60.4,-98}}));
        connect(HPO4.solution, solution.solution) annotation (Line(points={{20,-56},{22,-56},{22,-88},{60.4,-88},{60.4,-98}}));
        connect(PO4.solution, solution.solution) annotation (Line(points={{82,-72},{104,-72},{104,-108},{60.4,-108},{60.4,-98}}));
        connect(H.solution, solution.solution) annotation (Line(points={{10,-20},{10,-88},{60.4,-88},{60.4,-98}}));
      connect(H2O.solution, solution.solution) annotation (Line(
          points={{64,18},{104,18},{104,-108},{60.4,-108},{60.4,-98}},
          color={158,66,200}));
        connect(junctionN.outlet, H.inlet) annotation (Line(
            points={{-6,-8},{0,-8},{0,-10},{6,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction.products[1], H2PO4.inlet)
          annotation (Line(
            points={{-44,-48.25},{-44,-46},{-38,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction.products[2], junctionN.inlets[1])
          annotation (Line(
            points={{-44,-47.75},{-36,-47.75},{-36,-8.66667},{-26,-8.66667}},
            color={158,66,200},
            thickness=0.5));
        connect(H2PO4.outlet, chemicalReaction1.substrates[1]) annotation (Line(
            points={{-18,-46},{-18,-48},{-14,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction1.products[1], HPO4.inlet)
          annotation (Line(
            points={{6,-48.25},{12,-48.25},{12,-46},{16,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(HPO4.outlet, chemicalReaction2.substrates[1]) annotation (Line(
            points={{36,-46},{40,-46},{40,-48},{44,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction1.products[2], junctionN.inlets[2])
          annotation (Line(
            points={{6,-47.75},{10,-47.75},{10,-30},{-42,-30},{-42,-8},{-26,-8}},
            color={158,66,200},
            thickness=0.5));
        connect(H3PO4.outlet, chemicalReaction.substrates[1]) annotation (Line(
            points={{-72,-48},{-64,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction2.products[1], PO4.inlet)
          annotation (Line(
            points={{64,-48.25},{74,-48.25},{74,-62},{78,-62}},
            color={158,66,200},
            thickness=0.5));
        connect(chemicalReaction2.products[2], junctionN.inlets[3])
          annotation (Line(
            points={{64,-47.75},{72,-47.75},{72,-24},{-60,-24},{-60,-7.33333},{-26,-7.33333}},
            color={158,66,200},
            thickness=0.5));
        annotation ( Documentation(info="<html>
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=0.01, __Dymola_Algorithm="Dassl"));
      end Phosphate;

      model AcidBaseBufferTest
          extends Modelica.Icons.Example;

        Chemical.Boundaries.Buffer buffer(
          useInlet=true,
          substanceData(z=1.045),
          a_start=10^(-7.2),
          BufferValue=3) annotation (Placement(transformation(extent={{-50,4},{-30,24}})));

        Chemical.Solution simpleSolution annotation (Placement(transformation(extent={{-104,-100},{96,100}})));
        Chemical.Boundaries.ExternalMoleFraction externalMoleFraction(substanceData=Chemical.SubstancesOld.Proton_aqueous(), MoleFraction=10^(-7.1))
          annotation (Placement(transformation(extent={{0,-46},{20,-26}})));
        Chemical.Boundaries.Substance liquidWater(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{40,-80},{60,-60}})));
        Chemical.Processes.Diffusion process annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={16,20})));
      equation
        connect(buffer.solution, simpleSolution.solution) annotation (Line(
            points={{-46,4},{-26,4},{-26,-98},{56,-98}},
            color={127,127,0}));
        connect(liquidWater.solution, simpleSolution.solution)
          annotation (Line(points={{44,-80},{44,-98},{56,-98}}, color={127,127,0}));
        connect(externalMoleFraction.solution, simpleSolution.solution) annotation (Line(points={{4,-46},{4,-98},{56,-98}}, color={127,127,0}));
        connect(externalMoleFraction.outlet, process.inlet)
          annotation (Line(
            points={{20,-36},{26,-36},{26,-4},{16,-4},{16,10}},
            color={158,66,200},
            thickness=0.5));
        connect(process.outlet, buffer.inlet)
          annotation (Line(
            points={{16,30},{16,44},{-70,44},{-70,14},{-50,14}},
            color={158,66,200},
            thickness=0.5));
        annotation (                experiment(StopTime=0.05));
      end AcidBaseBufferTest;

      model CarbonDioxideInBlood
        import Chemical;
          extends Modelica.Icons.Example;

        parameter Real KC=1e-5;
        //e-6 "Slow down factor";
        Chemical.Solution blood_plasma(temperature_start=310.15) annotation (Placement(transformation(extent={{-100,4},{100,56}})));

        Chemical.Boundaries.Substance HCO3(
          use_mass_start=false,
          amountOfSubstance_start=0.024,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={18,24})));

        Chemical.Boundaries.ExternalIdealGas CO2_gas(
          substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(),
          PartialPressure(displayUnit="mmHg") = 5332.8954966,
          usePartialPressureInput=true) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,90})));
        Chemical.Processes.GasSolubility gasSolubility(
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0,
          substanceDataOut=Chemical.SubstancesOld.CarbonDioxide_aqueous()) annotation (Placement(transformation(extent={{-72,56},{-52,76}})));

        Chemical.Boundaries.Substance H2O(
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          use_mass_start=false,
          amountOfSubstance_start=51.6159) annotation (Placement(transformation(extent={{-28,14},{-48,34}})));
        Chemical.Boundaries.Substance Cl(
          substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.103) annotation (Placement(transformation(extent={{70,20},{50,40}})));

        Real  pH_p, pH_e;

        Modelica.Blocks.Sources.ContinuousClock clock(offset=5000)
          annotation (Placement(transformation(extent={{28,80},{8,100}})));
        Chemical.Boundaries.Substance others_P(
          substanceData=Chemical.Interfaces.Incompressible.SubstanceDataParameters(density=(1.024 - 0.933373)*1000/(1 - 0.936137), MolarWeight=(1.024 -
              0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024 - 0.0017)),
          use_mass_start=false,
          amountOfSubstance_start=0.1487) annotation (Placement(transformation(extent={{70,14},{90,34}})));
                  //References={"to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
        Chemical.Boundaries.Buffer H(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          BufferValue=0.0077,
          a_start=10^(-7.4)) "buffer value 7.7 mmol/L for plasma is from (O. Siggaard-Andersen 1995)"
          annotation (Placement(transformation(extent={{34,40},{52,58}})));

        Chemical.Solution blood_plasma1(ElectricGround=false, temperature_start=310.15)
                                                                 annotation (Placement(transformation(extent={{-96,-86},{104,-34}})));
        Chemical.Boundaries.Substance HCO3_E(
          use_mass_start=false,
          amountOfSubstance_start=0.0116,
          useInlet=true) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={12,-66})));

        Chemical.Boundaries.Substance CO2_E(
          use_mass_start=false,
          amountOfSubstance_start=0.0011,
          useInlet=true) "Free dissolved CO2 in red cells" annotation (Placement(transformation(extent={{-84,-82},{-64,-62}})));
        Chemical.Boundaries.Substance H2O_E(
          use_mass_start=false,
          amountOfSubstance_start=38.4008,
          useInlet=true) annotation (Placement(transformation(extent={{-54,-58},{-34,-38}})));
        Chemical.Boundaries.Substance Cl_E(
          use_mass_start=false,
          amountOfSubstance_start=0.0499,
          useInlet=true) annotation (Placement(transformation(extent={{52,-70},{72,-50}})));

        Chemical.Boundaries.Substance others_P1(
          substanceData=Chemical.Interfaces.Incompressible.SubstanceDataParameters(density=(1.024 - 0.933373)*1000/(1 - 0.936137), MolarWeight=(1.024 -
              0.933373)/(51.8*(1 - 0.994648) - 0.103 - 0.024 - 0.0017)),
          use_mass_start=false,
          amountOfSubstance_start=0.1444) annotation (Placement(transformation(extent={{74,-76},{94,-56}})));
              //References={"to reach plasma density 1024 kg/m3 and plasma volume 1 liter"},
        Chemical.Boundaries.Buffer H_E(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          BufferValue=0.063,
          a_start=10^(-7.2)) annotation (Placement(transformation(extent={{14,-52},{32,-34}})));

        Chemical.Processes.Membrane
                           Band3_Cl(
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={54,-18})));
        Chemical.Processes.Membrane
                           Band3_HCO3(
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-6,-20})));
        Chemical.Processes.Membrane Aquaporin(
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0)
                 annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-56,-16})));
        Chemical.Processes.Reaction HendersonHasselbalch1(
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          n_flow_0=0,
          productsSubstanceData={Chemical.SubstancesOld.Bicarbonate_blood(),Chemical.SubstancesOld.Proton_aqueous()},
          nS=2,
          nP=1) "K=10^(-6.103 + 3), dH=7.3 kJ/mol" annotation (Placement(transformation(extent={{-24,-58},{-4,-38}})));
        Chemical.Processes.Membrane ProtonExchanger(initN_flow=Chemical.Utilities.Types.InitializationMethods.state)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={36,-16})));
        Chemical.Topology.SplitterT1 splitterT1(redeclare package stateOfMatter = Chemical.Interfaces.IdealGas "Ideal Gas")
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-78,86})));
        Chemical.Boundaries.Substance CO2(
          use_mass_start=false,
          amountOfSubstance_start=0.00148,
          useInlet=true) "Free dissolved CO2 in plasma" annotation (Placement(transformation(extent={{-48,36},{-28,56}})));
        Chemical.Processes.Diffusion diffusion(initN_flow=Chemical.Utilities.Types.InitializationMethods.state)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-86,-18})));
        Chemical.Processes.FastGasSolubility gasSolubilitySubstance(substanceDataOut=Chemical.SubstancesOld.CarbonDioxide_aqueous(), redeclare package
            stateOfMatterIn = Chemical.Interfaces.IdealGas "Ideal Gas") annotation (Placement(transformation(extent={{-110,58},{-90,78}})));
      equation
        pH_p = -log10(H.substance.a);
        pH_e = -log10(H_E.substance.a);
      connect(H2O.solution, blood_plasma.solution)
        annotation (Line(points={{-32,14},{-32,0},{60,0},{60,4.52}},
                                                          color={127,127,0}));
      connect(Cl.solution, blood_plasma.solution) annotation (Line(
          points={{66,20},{66,12},{60,12},{60,4.52}},
          color={127,127,0}));
      connect(HCO3.solution, blood_plasma.solution) annotation (Line(
          points={{12,14},{12,12},{60,12},{60,8},{60,4},{60,4.52}},
          color={127,127,0}));
      connect(blood_plasma.solution, others_P.solution) annotation (Line(
          points={{60,4.52},{60,4},{60,8},{60,12},{74,12},{74,14}},
          color={127,127,0}));
      connect(clock.y, CO2_gas.partialPressure) annotation (Line(
          points={{7,90},{-10,90}},
          color={0,0,127}));
        connect(blood_plasma.solution, H.solution) annotation (Line(
            points={{60,4.52},{60,12},{37.6,12},{37.6,40}},
            color={127,127,0}));
        connect(CO2_gas.solution, blood_plasma.solution) annotation (Line(points={{-14,100},{-40,100},{-40,108},{104,108},{104,0},{60,0},{60,4.52}},
                                                                                                                                   color={127,127,0}));
        connect(CO2_E.solution, blood_plasma1.solution) annotation (Line(points={{-80,-82},{-80,-100},{64,-100},{64,-85.48}},
                                                                                                                            color={127,127,0}));
        connect(H2O_E.solution, blood_plasma1.solution) annotation (Line(points={{-50,-58},{-50,-100},{64,-100},{64,-85.48}},
                                                                                                                  color={127,127,0}));
        connect(Cl_E.solution, blood_plasma1.solution) annotation (Line(points={{56,-70},{56,-78},{64,-78},{64,-85.48}}, color={127,127,0}));
        connect(HCO3_E.solution, blood_plasma1.solution) annotation (Line(points={{6,-76},{6,-90},{64,-90},{64,-85.48}},   color={127,127,0}));
        connect(blood_plasma1.solution, others_P1.solution) annotation (Line(points={{64,-85.48},{64,-78},{78,-78},{78,-76}}, color={127,127,0}));
        connect(blood_plasma1.solution, H_E.solution) annotation (Line(points={{64,-85.48},{64,-78},{38,-78},{38,-52},{17.6,-52}}, color={127,127,0}));
        connect(Band3_Cl.outlet, Cl_E.inlet) annotation (Line(
            points={{54,-28},{54,-44},{52,-44},{52,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(Band3_Cl.inlet, Cl.outlet) annotation (Line(
            points={{54,-8},{54,30},{50,30}},
            color={158,66,200},
            thickness=0.5));
        connect(Band3_HCO3.outlet, HCO3.inlet) annotation (Line(
            points={{-6,-10},{-6,24},{8,24}},
            color={158,66,200},
            thickness=0.5));
        connect(H2O.outlet, Aquaporin.inlet) annotation (Line(
            points={{-48,24},{-48,-6},{-56,-6}},
            color={158,66,200},
            thickness=0.5));
        connect(Aquaporin.outlet, H2O_E.inlet) annotation (Line(
            points={{-56,-26},{-56,-38},{-54,-38},{-54,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(H2O_E.outlet, HendersonHasselbalch1.substrates[1])
          annotation (Line(
            points={{-34,-48},{-34,-48.25},{-24,-48.25}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2_E.outlet, HendersonHasselbalch1.substrates[2])
          annotation (Line(
            points={{-64,-72},{-48,-72},{-48,-70},{-24,-70},{-24,-47.75}},
            color={158,66,200},
            thickness=0.5));
        connect(HCO3_E.outlet, Band3_HCO3.inlet)
          annotation (Line(
            points={{22,-66},{32,-66},{32,-54},{4,-54},{4,-40},{-6,-40},{-6,-30}},
            color={158,66,200},
            thickness=0.5));
        connect(H_E.outlet, ProtonExchanger.inlet) annotation (Line(
            points={{32,-43},{32,-42},{36,-42},{36,-26}},
            color={158,66,200},
            thickness=0.5));
        connect(ProtonExchanger.outlet, H.inlet) annotation (Line(
            points={{36,-6},{34,-6},{34,49}},
            color={158,66,200},
            thickness=0.5));
        connect(diffusion.outlet, CO2_E.inlet) annotation (Line(
            points={{-86,-28},{-86,-72},{-84,-72}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.solution, blood_plasma.solution) annotation (Line(points={{-44,36},{-44,64},{104,64},{104,0},{60,0},{60,4.52}},
                                                                                                              color={127,127,0}));
        connect(gasSolubility.outlet, CO2.inlet) annotation (Line(
            points={{-52,66},{-48,66},{-48,46}},
            color={158,66,200},
            thickness=0.5));
        connect(gasSolubilitySubstance.outlet, diffusion.inlet) annotation (Line(
            points={{-90,68},{-86,68},{-86,-8}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2_gas.outlet, splitterT1.inlet)
          annotation (Line(
            points={{-30,90},{-58,90},{-58,106},{-78,106},{-78,96}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletA, gasSolubility.inlet) annotation (Line(
            points={{-68,86},{-68,66},{-72,66}},
            color={158,66,200},
            thickness=0.5));
        connect(gasSolubilitySubstance.inlet, splitterT1.outletB)
          annotation (Line(
            points={{-110,68},{-118,68},{-118,86},{-88,86}},
            color={158,66,200},
            thickness=0.5));
        connect(gasSolubilitySubstance.solution, blood_plasma.solution) annotation (Line(points={{-89.8,63.8},{-120,63.8},{-120,108},{104,108},{104,0},{60,0},{60,
                4.52}},                                                                                                                     color={127,127,0}));
        connect(blood_plasma1.solution, diffusion.solution)
          annotation (Line(points={{64,-85.48},{-48,-85.48},{-48,-86},{-90.2,-86},{-90.2,-28.2}}, color={127,127,0}));
        connect(Aquaporin.solution, blood_plasma1.solution) annotation (Line(points={{-60.2,-26.2},{-60.2,-85.48},{64,-85.48}}, color={127,127,0}));
        connect(HendersonHasselbalch1.products[1], HCO3_E.inlet)
          annotation (Line(
            points={{-4,-48},{0,-48},{0,-66},{2,-66}},
            color={158,66,200},
            thickness=0.5));
        connect(Band3_HCO3.solution, blood_plasma.solution) annotation (Line(points={{-10.2,-9.8},{-10.2,4.52},{60,4.52}}, color={127,127,0}));
        connect(Band3_Cl.solution, blood_plasma1.solution) annotation (Line(points={{49.8,-28.2},{52,-100},{64,-100},{64,-85.48}}, color={127,127,0}));
        annotation ( Documentation(info="<html>
<p>The mature red blood cell (erythrocyte) is the simplest cell in the human body. Its primary function is the transportation of blood gases, such as oxygen O<sub>2</sub> (from the lungs to tissues) and carbon dioxide CO<sub>2</sub> (from tissues to the lungs). The chemical processes behind the gases&rsquo; transportation are complex because the capacity of water to transport their freely dissolved forms is very low. To transport sufficient amounts of O<sub>2</sub> and CO<sub>2</sub>, the gases must be chemically bound to hemoglobin such as described in (Matej&aacute;k, et al., 2015) and/or transported as different substances, which can be present in water in much higher concentrations than their freely dissolved forms allow. Therefore, to transport a sufficient amount of CO<sub>2</sub>, it must be changed to HCO<sub>3</sub><sup>-</sup> using the chemical reaction: </p>
<table width=100%><tr>
<th><p align=\"center\"><b>CO<sub>2</sub> + H<sub>2</sub>O &lt;-&gt; HCO<sub>3</sub><sup>-</sup> + H<sup>+</sup></b></p></th>
<td><p>(1)</p></td>
</tr>
</table>
<p><br>This reaction takes place mainly inside the red blood cell, because only here it is presented with the enzyme carbonic anhydrase. Therefore, the increase of total carbon dioxide content of blood in tissues and its decrease in lungs are always connected with the chloride shift between blood plasma and the intracellular fluid of erythrocytes, as represented in followin Figure: </p>
<p><img src=\"modelica://Chemical/Resources/Images/Examples/CO2inBlood.png\"/></p>
<p>Figure: Chloride shift with carbon dioxide hydration with assumption of non-bicarbonate linear acid-base buffering properties of plasma and erythrocytes. </p>
<p><br>The blood plasma and intracellular fluid are divided by the cellular membrane composed of a special, very compact lipid double-layer. A lipophobic compound (not soluble in lipids) cannot cross the membrane without special proteins called membrane channels. Even water molecules must have membrane channels (called aquaporins) in order to cross the cellular membrane. In addition, the chloride shift (also known as the Hamburger shift) is exchanging an aqueous chloride Cl<sup>-</sup> for an aqueous bicarbonate HCO<sub>3</sub><sup>-</sup> in both directions across the cellular membranes of red blood cells using the membrane channel &ldquo;Band 3&rdquo;. Each passive membrane channel only allows the equilibration of the electrochemical potentials of the specific permeable ions on both sides of membrane. The different electric potentials on each side of membrane allow their different concentrations to achieve equilibrium. </p>
<p>Conversely, the solution&rsquo;s equilibrium of different ions&rsquo; compositions on both sides of the membrane creates the measurable electric membrane potential. This process is not so intuitive, because even though neither solution needs to have an electric charge, there can be a non-zero electric potential for permeable ions. This potential for permeable ions at equilibrium is called the Nernst membrane potential and, in the Chemical library, it is a direct mathematical result of the equality of the electrochemical potential of the ion in both solutions. </p>
<p>The intracellular solution must be set at the possible nonzero electric potential (ElectricalGround=false) because, as a result, the membrane potential of the erythrocytes is calculated as -12mV, which agrees with experimental data by Gedde and Huestis (Gedde and Huestis, 1997) in the electrolytes&rsquo; setting by Raftos et al. (Raftos, et al., 1990). </p>
<p>In this way, it is possible to model more complex processes of a membrane where chemical reactions of active membrane channels or membrane receptors can both be used.&nbsp; </p>
<p><br>CO2 in blood with linear H+ non-bicarbonates buffering without binding to hemoglobin.</p>
<p>The buffer values 0.063 mmol/L commes from Siggaard-Andersen.</p>
</html>",      revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1500, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          __Dymola_experimentSetupOutput,
          __Dymola_experimentFlags(
            Advanced(
              EvaluateAlsoTop=false,
              GenerateAnalyticJacobian=false,
              OutputModelicaCode=false),
            Evaluate=true,
            OutputCPUtime=false,
            OutputFlatModelica=false));
      end CarbonDioxideInBlood;

      model RedCellMembrane
        import Chemical;
        import Chemical;
       // import Chemical;
        extends Modelica.Icons.Example;

        parameter Real KC=1e-5;
        //e-6 "Slow down factor";
        Chemical.Solution blood_erythrocytes(ElectricGround=false) annotation (Placement(transformation(extent={{-180,-100},{180,-10}})));
        Chemical.Solution blood_plasma annotation (Placement(transformation(extent={{-180,12},{180,100}})));

        Chemical.Boundaries.Substance HCO3(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.024) annotation (Placement(transformation(extent={{10,-10},{-10,10}}, origin={-18,30})));

        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=51.8*0.994648/55.508)
          annotation (Placement(transformation(extent={{-146,44},{-166,64}})));
        Chemical.Boundaries.Substance HCO3_E(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Bicarbonate_blood(),
          use_mass_start=false,
          amountOfSubstance_start=0.0116) annotation (Placement(transformation(extent={{-28,-38},{-8,-18}})));
        Chemical.Boundaries.Substance H2O_E(
          use_mass_start=true,
          mass_start=0.68,
          useInlet=true) annotation (Placement(transformation(extent={{-164,-38},{-144,-18}})));
        Chemical.Boundaries.Substance Cl_E(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.0499) annotation (Placement(transformation(extent={{22,-36},{42,-16}})));
        Chemical.Boundaries.Substance Cl(
          substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.103) annotation (Placement(transformation(extent={{-4,20},{16,40}})));
        Chemical.Boundaries.Substance albumin(
          substanceData(
            MolarWeight=66.463,
            z=-17,
            density=1080),
          use_mass_start=false,
          amountOfSubstance_start=0.0007) annotation (Placement(transformation(extent={{112,76},{92,96}})));
        Real pH_p; //pH_e,

        Chemical.Processes.Membrane Aquapirin
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-168,0})));
        Chemical.Processes.Membrane Band3 annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-6,-2})));
        Chemical.Processes.Membrane Band3_ annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={18,0})));
        Chemical.Boundaries.Substance permeableUncharged(use_mass_start=false, amountOfSubstance_start=0.0118)
          annotation (Placement(transformation(extent={{166,20},{146,40}})));
        Chemical.Boundaries.Substance permeableUncharged_E(
          useInlet=true,
          substanceData(MolarWeight=0.1),
          use_mass_start=false,
          amountOfSubstance_start=0.00903) annotation (Placement(transformation(extent={{144,-38},{164,-18}})));
        Chemical.Boundaries.Substance chargedImpermeable_E(
          substanceData(MolarWeight=1),
          use_mass_start=false,
          amountOfSubstance_start=0.0165) annotation (Placement(transformation(extent={{144,-62},{164,-42}})));
        Chemical.Processes.Membrane leak annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={140,0})));
        Chemical.Boundaries.Substance Lac_E(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00062) annotation (Placement(transformation(extent={{76,-38},{56,-18}})));
        Chemical.Boundaries.Substance Lac(
          substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00131) annotation (Placement(transformation(extent={{56,20},{76,40}})));
        Chemical.Processes.Membrane MCT_ "Monocarboxylate transporters"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={80,0})));
        Chemical.Boundaries.Substance H(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=10^(-7.4)) "H+ in plasma" annotation (Placement(transformation(extent={{50,20},{30,40}})));
        Chemical.Processes.Membrane MCT "Monocarboxylate transporters"
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={52,0})));
        Chemical.Boundaries.Substance CO2(
          substanceData=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00167) "free dissolved unbound CO2" annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        Chemical.Boundaries.Substance CO2_E(
          use_mass_start=false,
          amountOfSubstance_start=0.00125,
          useInlet=true) "free dissolved unbound CO2" annotation (Placement(transformation(extent={{-38,-38},{-58,-18}})));
        Chemical.Processes.Membrane freeCO2 annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-38,2})));
        Chemical.Boundaries.Substance O2(
          substanceData=Chemical.SubstancesOld.Oxygen_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000167) "free dissolved undound oxygen" annotation (Placement(transformation(extent={{96,20},{116,40}})));
        Chemical.Processes.Membrane freeO2 annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={118,0})));
        Chemical.Boundaries.Substance O2_E(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Oxygen_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000125) "free dissolved undound O2" annotation (Placement(transformation(extent={{116,-38},{96,-18}})));
        Chemical.Boundaries.Substance K(
          substanceData=Chemical.SubstancesOld.Potassium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.004) annotation (Placement(transformation(extent={{-100,20},{-120,40}})));
        Chemical.Boundaries.Substance Na(
          substanceData=Chemical.SubstancesOld.Sodium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.138) annotation (Placement(transformation(extent={{-124,20},{-144,40}})));
        Chemical.Boundaries.Substance Na_E(
          substanceData=Chemical.SubstancesOld.Sodium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.007) annotation (Placement(transformation(extent={{-118,-38},{-138,-18}})));
        Chemical.Boundaries.Substance K_E(
          substanceData=Chemical.SubstancesOld.Potassium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.096) annotation (Placement(transformation(extent={{-112,-38},{-92,-18}})));
        Chemical.Boundaries.Substance H2PO4_E(
          substanceData=Chemical.SubstancesOld.DihydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000175) annotation (Placement(transformation(extent={{-84,-38},{-64,-18}})));
        Chemical.Boundaries.Substance ADP_E(
          substanceData(z=-3),
          use_mass_start=false,
          amountOfSubstance_start=9.6e-05) annotation (Placement(transformation(extent={{-114,-62},{-94,-42}})));
        Chemical.Boundaries.Substance ATP_E(
          substanceData(
            z=-4,
            DfH=16700,
            DfG=30500),
          use_mass_start=false,
          amountOfSubstance_start=0.00128) annotation (Placement(transformation(extent={{-146,-62},{-166,-42}})));
            //References={"http://www.wiley.com/college/pratt/0471393878/student/review/thermodynamics/7_relationship.html"}

        Chemical.Boundaries.Substance HPO4_E(
          substanceData=Chemical.SubstancesOld.HydrogenPhosphate_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.000495) annotation (Placement(transformation(extent={{-84,-62},{-64,-42}})));
        Chemical.Boundaries.Substance globulins(
          substanceData(
            MolarWeight=34,
            z=-2.43,
            density=1080),
          use_mass_start=false,
          amountOfSubstance_start=0.00082) annotation (Placement(transformation(extent={{150,76},{130,96}})));
        Chemical.Boundaries.Substance Ca(
          substanceData=Chemical.SubstancesOld.Calcium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00175) "Ca2+" annotation (Placement(transformation(extent={{-78,20},{-98,40}})));
        Chemical.Boundaries.Substance Mg(
          substanceData=Chemical.SubstancesOld.Magnesium_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=0.00108) "Mg2+" annotation (Placement(transformation(extent={{-112,-84},{-92,-64}})));
        Chemical.Boundaries.Substance DPG(
          substanceData(
            MolarWeight=0.266,
            z=-2.2,
            density=1000),
          use_mass_start=false,
          amountOfSubstance_start=0.0051) annotation (Placement(transformation(extent={{128,-94},{108,-74}})));
        Chemical.Boundaries.Substance GSH(
          substanceData(
            MolarWeight=0.2,
            z=-1,
            density=1000),
          use_mass_start=false,
          amountOfSubstance_start=0.00223) annotation (Placement(transformation(extent={{164,-94},{144,-74}})));
        Chemical.Processes.Reaction HendersonHasselbalch(nS=2, nP=2) "K=10^(-6.103 + 3), dH=7.3 kJ/mol"
          annotation (Placement(transformation(extent={{-24,-58},{-4,-78}})));

        Chemical.Boundaries.Buffer Hemoglobin(
          useInlet=true,
          substanceData(z=1.045),
          a_start=10^(-7.2),
          BufferValue=3) annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={34,-68})));
      equation
        pH_p = -log10(H.substance.a);
      //  pH_e = -log10(H_E.a);
      connect(H2O.solution, blood_plasma.solution)
        annotation (Line(points={{-150,44},{-150,44},{-150,26},{-150,26},{-150,20},{108,
                20},{108,12.88}}, color={127,127,0}));
      connect(Cl.solution, blood_plasma.solution) annotation (Line(
          points={{0,20},{0,16},{0,12.88},{108,12.88}},
          color={127,127,0}));
        connect(H2O_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-160,-38},{108,-38},{108,-99.1}},
                                                      color={127,127,0}));
        connect(Cl_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{26,-36},{28,-36},{28,-99.1},{108,-99.1}},
            color={127,127,0}));
        connect(HCO3_E.solution, blood_erythrocytes.solution) annotation (Line(
              points={{-24,-38},{108,-38},{108,-99.1}},
                                                    color={127,127,0}));
      connect(HCO3.solution, blood_plasma.solution) annotation (Line(
          points={{-12,20},{108,20},{108,12.88}},
          color={127,127,0}));
        connect(blood_plasma.solution, permeableUncharged.solution) annotation (Line(
            points={{108,12.88},{108,20},{162,20}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution, permeableUncharged_E.solution)
          annotation (Line(
            points={{108,-99.1},{108,-38},{148,-38}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution,chargedImpermeable_E. solution)
          annotation (Line(
            points={{108,-99.1},{108,-38},{140,-38},{140,-62},{148,-62}},
            color={127,127,0}));
        connect(Lac.solution, blood_plasma.solution) annotation (Line(
            points={{60,20},{108,20},{108,12.88}},
            color={127,127,0}));
        connect(blood_erythrocytes.solution, Lac_E.solution) annotation (Line(
            points={{108,-99.1},{108,-38},{72,-38}},
            color={127,127,0}));
        connect(blood_plasma.solution, H.solution) annotation (Line(
            points={{108,12.88},{108,20},{46,20}},
            color={127,127,0}));
        connect(blood_plasma.solution, CO2.solution) annotation (Line(
            points={{108,12.88},{108,20},{-56,20}},
            color={127,127,0}));
        connect(CO2_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-42,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(blood_plasma.solution, O2.solution) annotation (Line(
            points={{108,12.88},{108,20},{100,20}},
            color={127,127,0}));
        connect(O2_E.solution, blood_erythrocytes.solution) annotation (Line(
            points={{112,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(H2O.solution, K.solution) annotation (Line(
            points={{-150,44},{-150,20},{-104,20}},
            color={158,66,200}));
        connect(H2O.solution, Na.solution) annotation (Line(
            points={{-150,44},{-150,20},{-128,20}},
            color={158,66,200}));
        connect(H2O_E.solution, Na_E.solution) annotation (Line(
            points={{-160,-38},{-122,-38}},
            color={158,66,200}));
        connect(H2O_E.solution, K_E.solution) annotation (Line(
            points={{-160,-38},{-108,-38}},
            color={158,66,200}));
        connect(H2O_E.solution, H2PO4_E.solution) annotation (Line(
            points={{-160,-38},{-80,-38}},
            color={127,127,0}));
        connect(ADP_E.solution, K_E.solution) annotation (Line(
            points={{-110,-62},{-110,-38},{-108,-38}},
            color={158,66,200}));
        connect(ATP_E.solution, Na_E.solution) annotation (Line(
            points={{-150,-62},{-122,-62},{-122,-38}},
            color={127,127,0}));
        connect(HPO4_E.solution, H2PO4_E.solution) annotation (Line(
            points={{-80,-62},{-110,-62},{-110,-38},{-80,-38}},
            color={127,127,0}));
        connect(Ca.solution, CO2.solution) annotation (Line(
            points={{-82,20},{-82,20},{-56,20}},
            color={127,127,0}));
        connect(Mg.solution, blood_erythrocytes.solution) annotation (Line(
            points={{-108,-84},{-108,-38},{108,-38},{108,-99.1}},
            color={127,127,0}));
        connect(DPG.solution, permeableUncharged_E.solution) annotation (Line(
            points={{124,-94},{140,-94},{140,-38},{148,-38}},
            color={127,127,0}));
        connect(GSH.solution, permeableUncharged_E.solution) annotation (Line(
            points={{160,-94},{140,-94},{140,-38},{148,-38}},
            color={127,127,0}));
        connect(Hemoglobin.solution, blood_erythrocytes.solution) annotation (Line(
              points={{28,-78},{28,-88},{108,-88},{108,-99.1}}, color={127,127,
                0}));

        connect(albumin.solution, blood_plasma.solution) annotation (Line(
            points={{108,76},{126,76},{126,20},{108,20},{108,12.88}},
            color={127,127,0},
            smooth=Smooth.None));
        connect(globulins.solution, blood_plasma.solution) annotation (Line(
            points={{146,76},{126,76},{126,20},{108,20},{108,12.88}},
            color={127,127,0},
            smooth=Smooth.None));
        connect(H2O.outlet, Aquapirin.inlet) annotation (Line(
            points={{-166,54},{-168,54},{-168,10}},
            color={158,66,200},
            thickness=0.5));
        connect(Aquapirin.outlet, H2O_E.inlet) annotation (Line(
            points={{-168,-10},{-168,-28},{-164,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.outlet, freeCO2.inlet) annotation (Line(
            points={{-40,30},{-38,30},{-38,12}},
            color={158,66,200},
            thickness=0.5));
        connect(freeCO2.outlet, CO2_E.inlet) annotation (Line(
            points={{-38,-8},{-38,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(H2O_E.outlet, HendersonHasselbalch.substrates[1])
          annotation (Line(
            points={{-144,-28},{-142,-28},{-142,-90},{-48,-90},{-48,-67.75},{-24,-67.75}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2_E.outlet, HendersonHasselbalch.substrates[2])
          annotation (Line(
            points={{-58,-28},{-60,-28},{-60,-68.25},{-24,-68.25}},
            color={158,66,200},
            thickness=0.5));
        connect(Band3.outlet, HCO3.inlet) annotation (Line(
            points={{-6,8},{-6,30},{-8,30}},
            color={158,66,200},
            thickness=0.5));
        connect(HendersonHasselbalch.products[1], HCO3_E.inlet)
          annotation (Line(
            points={{-4,-67.75},{0,-67.75},{0,-54},{-28,-54},{-28,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(HendersonHasselbalch.products[2], Hemoglobin.inlet)
          annotation (Line(
            points={{-4,-68.25},{-4,-70},{24,-70},{24,-68}},
            color={158,66,200},
            thickness=0.5));
        connect(Hemoglobin.outlet, MCT.inlet) annotation (Line(
            points={{44,-68},{52,-68},{52,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(MCT.outlet, H.inlet) annotation (Line(
            points={{52,10},{52,30},{50,30}},
            color={158,66,200},
            thickness=0.5));
        connect(Band3_.outlet, Cl_E.inlet) annotation (Line(
            points={{18,-10},{18,-26},{22,-26}},
            color={158,66,200},
            thickness=0.5));
        connect(Cl.outlet, Band3_.inlet) annotation (Line(
            points={{16,30},{18,30},{18,10}},
            color={158,66,200},
            thickness=0.5));
        connect(Lac.outlet, MCT_.inlet) annotation (Line(
            points={{76,30},{80,30},{80,10}},
            color={158,66,200},
            thickness=0.5));
        connect(MCT_.outlet, Lac_E.inlet) annotation (Line(
            points={{80,-10},{80,-28},{76,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(O2.outlet, freeO2.inlet) annotation (Line(
            points={{116,30},{118,30},{118,10}},
            color={158,66,200},
            thickness=0.5));
        connect(freeO2.outlet, O2_E.inlet) annotation (Line(
            points={{118,-10},{120,-10},{120,-28},{116,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(permeableUncharged.outlet, leak.inlet) annotation (Line(
            points={{146,30},{140,30},{140,10}},
            color={158,66,200},
            thickness=0.5));
        connect(leak.outlet, permeableUncharged_E.inlet) annotation (Line(
            points={{140,-10},{140,-28},{144,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(HCO3_E.outlet, Band3.inlet) annotation (Line(
            points={{-8,-28},{-6,-28},{-6,-12}},
            color={158,66,200},
            thickness=0.5));
        connect(Aquapirin.solution, blood_erythrocytes.solution) annotation (Line(points={{-172.2,-10.2},{-172.2,-99.1},{108,-99.1}}, color={127,127,0}));
        connect(freeCO2.solution, blood_erythrocytes.solution) annotation (Line(points={{-42.2,-8.2},{-42.2,-99.1},{108,-99.1}}, color={127,127,0}));
        annotation ( Documentation(info="<html>
<p>Blood eqiulibrium across erythrocyte membrane bewteen blood plasma and intracellular fluid of erythrocytes.</p>
<p>Data of blood status are from:</p>
<p>Raftos, J.E., Bulliman, B.T. and Kuchel, P.W. Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. <i>The Journal of general physiology</i> 1990;95(6):1183-1204. </p>
</html>",      revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),experiment(StopTime=1e-08, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-100},
                  {180,100}})),
          __Dymola_experimentSetupOutput,
          __Dymola_experimentFlags(
            Advanced(
              EvaluateAlsoTop=false,
              GenerateAnalyticJacobian=false,
              OutputModelicaCode=false),
            Evaluate=true,
            OutputCPUtime=false,
            OutputFlatModelica=false));
      end RedCellMembrane;

      package Dev

        model NaKATPase
          Solution ICF(ElectricGround=false) annotation (Placement(transformation(extent={{-100,-100},{-20,98}})));
          Solution ECF annotation (Placement(transformation(extent={{28,-100},{100,100}})));
          Chemical.Boundaries.Substance water_ICF(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
            annotation (Placement(transformation(extent={{-82,70},{-62,90}})));
          Chemical.Boundaries.Substance water_EFC(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
            annotation (Placement(transformation(extent={{58,56},{38,76}})));
          Chemical.Processes.Reaction NaKATPase(
            nS=4,
            nP=4,
            s={3,2,1,1},
            p={3,2,1,1}) annotation (Placement(transformation(
                extent={{-10,10},{10,-10}},
                rotation=270,
                origin={0,8})));
          Chemical.Boundaries.Substance Na_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.01,
            substanceData=Chemical.SubstancesOld.Sodium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-72,18},{-52,38}})));
          Chemical.Boundaries.Substance K_ICF(
            useInlet=true,
            amountOfSubstance_start(displayUnit="mmol") = 0.14,
            substanceData=Chemical.SubstancesOld.Potassium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-56,-18},{-76,2}})));
          Chemical.Boundaries.Substance ATP4_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.0038,
            substanceData=Chemical.SubstancesOld.ATP4_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-72,42},{-52,62}})));
          Chemical.Boundaries.Substance Na_ECF(
            useInlet=true,
            amountOfSubstance_start(displayUnit="mmol") = 0.14,
            substanceData=Chemical.SubstancesOld.Sodium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{42,-18},{62,2}})));
          Chemical.Boundaries.Substance K_ECF(
            amountOfSubstance_start(displayUnit="mmol") = 0.005,
            substanceData=Chemical.SubstancesOld.Potassium_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{58,28},{38,48}})));
          Chemical.Boundaries.Substance ADP3(
            useInlet=true,
            substanceData=Chemical.SubstancesOld.ADP3_aqueous(),
            use_mass_start=false,
            amountOfSubstance_start=0.005/150) annotation (Placement(transformation(extent={{-56,-42},{-76,-22}})));
          Chemical.Boundaries.Substance H2PO4_ICF(
            useInlet=true,
            amountOfSubstance_start(displayUnit="mmol") = 0.036,
            substanceData=Chemical.SubstancesOld.DihydrogenPhosphate_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-56,-68},{-76,-48}})));
          inner Modelica.Fluid.System system(T_ambient=310.15)
            annotation (Placement(transformation(extent={{-6,-90},{14,-70}})));
          Chemical.Boundaries.Substance Cl_ICF(
            amountOfSubstance_start(displayUnit="mmol") = 0.0987,
            substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{-82,-92},{-62,-72}})));
          Chemical.Boundaries.Substance Cl_ECF(
            amountOfSubstance_start(displayUnit="mmol") = 0.145,
            substanceData=Chemical.SubstancesOld.Chloride_aqueous(),
            use_mass_start=false) annotation (Placement(transformation(extent={{50,-72},{70,-52}})));
        equation
          connect(water_ICF.solution, ICF.solution) annotation (Line(points={{-78,70},
                  {-78,-100},{-36,-100},{-36,-98.02}}, color={127,127,0}));
          connect(water_EFC.solution, ECF.solution) annotation (Line(points={{54,56},
                  {54,-98},{85.6,-98}}, color={127,127,0}));
          connect(water_EFC.solution, Na_ECF.solution)
            annotation (Line(points={{54,56},{54,-18},{46,-18}},
                                                        color={127,127,0}));
          connect(water_EFC.solution, K_ECF.solution)
            annotation (Line(points={{54,56},{54,28}}, color={127,127,0}));
          connect(water_ICF.solution, Na_ICF.solution)
            annotation (Line(points={{-78,70},{-78,18},{-68,18}},
                                                         color={127,127,0}));
          connect(water_ICF.solution, K_ICF.solution)
            annotation (Line(points={{-78,70},{-78,-18},{-60,-18}},
                                                         color={127,127,0}));
          connect(water_ICF.solution, ADP3.solution)
            annotation (Line(points={{-78,70},{-78,-42},{-60,-42}},
                                                          color={127,127,0}));
          connect(water_ICF.solution, H2PO4_ICF.solution)
            annotation (Line(points={{-78,70},{-78,-68},{-60,-68}},
                                                          color={127,127,0}));
          connect(Na_ICF.outlet, NaKATPase.substrates[1])
            annotation (Line(points={{-52,28},{0.375,28},{0.375,18}},
                                                                color={158,66,200}));
          connect(K_ECF.outlet, NaKATPase.substrates[2])
            annotation (Line(points={{38,38},{0.125,38},{0.125,18}},
                                                               color={158,66,200}));
          connect(NaKATPase.products[1], Na_ECF.inlet) annotation (Line(points={{0.375,-2},{-2,-2},{-2,-8},{42,-8}},
                                                  color={158,66,200}));
          connect(NaKATPase.products[2], K_ICF.inlet) annotation (Line(points={{0.125,-2},{0,-2},{0,-8},{-56,-8}},
                                                             color={158,66,200}));
          connect(ATP4_ICF.outlet, NaKATPase.substrates[3]) annotation (Line(points={{-52,52},{-0.125,52},{-0.125,18}},
                                                             color={158,66,200}));
          connect(water_ICF.outlet, NaKATPase.substrates[4])
            annotation (Line(points={{-62,80},{-0.375,80},{-0.375,18}},
                                                              color={158,66,200}));
          connect(NaKATPase.products[3], ADP3.inlet) annotation (Line(points={{-0.125,-2},{2,-2},{2,-32},{-56,-32}},
                                             color={158,66,200}));
          connect(H2PO4_ICF.inlet, NaKATPase.products[4])
            annotation (Line(points={{-56,-58},{-0.375,-58},{-0.375,-2}},
                                                                color={158,66,200}));
          connect(water_ICF.solution, Cl_ICF.solution)
            annotation (Line(points={{-78,70},{-78,-92}}, color={127,127,0}));
          connect(water_EFC.solution, Cl_ECF.solution)
            annotation (Line(points={{54,56},{54,-72}}, color={127,127,0}));
          connect(ATP4_ICF.solution, ICF.solution) annotation (Line(points={{-68,42},{-104,42},{-104,-104},{-36,-104},{-36,-98.02}}, color={127,127,0}));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end NaKATPase;
      end Dev;

      model AlbuminTitration "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

        Chemical.Solution solution(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible)
          annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

        constant Integer n=218 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]=cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
          "acid dissociation constants";
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Boundaries.Substance A[n](
          substanceData(each z=-1),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "deprotonated acid groups" annotation (Placement(transformation(extent={{26,-16},{6,4}})));
        Chemical.Processes.Reaction react[n](
          initN_flow=Chemical.Utilities.Types.InitializationMethods.state,
          each nS=2,
          each nP=1) annotation (Placement(transformation(extent={{-24,-2},{-44,18}})));

        Chemical.Boundaries.Substance HA[n](
          useInlet=true,
          substanceData(DfG=DfG),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "protonated acid groups" annotation (Placement(transformation(extent={{-58,-2},{-78,18}})));

        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,-68})));
        Chemical.Boundaries.ExternalMoleFraction H(substanceData=Chemical.SubstancesOld.Proton_aqueous(), MoleFraction=10^(-7.4))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,42})));
        inner DropOfCommons dropOfCommons(L=1e-5) annotation (Placement(transformation(extent={{62,68},{82,88}})));
      equation

        for i in 1:n loop
          connect(HA[i].solution, solution.solution) annotation (Line(
            points={{-62,-2},{-62,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(A[i].solution, solution.solution) annotation (Line(
            points={{22,-16},{22,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(react[i].substrates[1], A[i].outlet) annotation (Line(
            points={{-24,7.75},{-12,7.75},{-12,-6},{6,-6}},
            color={107,45,134},
            thickness=1));
          connect(H.outlet, react[i].substrates[2]) annotation (Line(points={{10,42},{-8,42},{-8,8.25},{-24,8.25}},
                                color={158,66,200}));
          connect(HA[i].inlet, react[i].products[1]) annotation (Line(
            points={{-58,8},{-50,8},{-50,8},{-44,8}},
            color={107,45,134},
            thickness=1));
        end for;

        connect(solution.solution, H2O.solution) annotation (Line(
          points={{56,-98},{56,-78}},
          color={127,127,0}));

        connect(H.solution, solution.solution) annotation (Line(points={{26,52},{26,78},{-94,78},{-94,-86},{56,-86},{56,-98}}, color={127,127,0}));
        annotation ( Documentation(revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"));
      end AlbuminTitration;

      model AlbuminTitration1 "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

        Chemical.Solution solution(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible)
          annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

        constant Integer n=1 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]={8.5}; //cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
          //"acid dissociation constants";
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Boundaries.Substance A[n](
          useInlet=true,
          substanceData(each z=-1),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "deprotonated acid groups" annotation (Placement(transformation(extent={{6,-16},{26,4}})));
        Chemical.Processes.Reaction react[n](
          each nS=1,
          each nP=2) annotation (Placement(transformation(extent={{-44,-2},{-24,18}})));

        Chemical.Boundaries.Substance HA[n](
          substanceData(DfG=DfG),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "protonated acid groups" annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));

        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,-68})));
        Chemical.Boundaries.ExternalMoleFraction H(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          MoleFraction=10^(-7.4)) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={20,42})));
      equation

        for i in 1:n loop
          connect(HA[i].solution, solution.solution) annotation (Line(
            points={{-74,-2},{-74,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(A[i].solution, solution.solution) annotation (Line(
            points={{10,-16},{10,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(react[i].products[1], A[i].inlet) annotation (Line(
            points={{-24,7.75},{-12,7.75},{-12,-6},{6,-6}},
            color={107,45,134},
            thickness=1));
          connect(H.inlet, react[i].products[2]) annotation (Line(points={{10,42},{-8,42},{-8,8.25},{-24,8.25}},
                                color={158,66,200}));
          connect(HA[i].outlet, react[i].substrates[1]) annotation (Line(
            points={{-58,8},{-50,8},{-50,8},{-44,8}},
            color={107,45,134},
            thickness=1));
        end for;

        connect(solution.solution, H2O.solution) annotation (Line(
          points={{56,-98},{56,-78}},
          color={127,127,0}));

        connect(H.solution, solution.solution) annotation (Line(points={{14,52},{14,78},{-94,78},{-94,-86},{56,-86},{56,-98}}, color={127,127,0}));
        annotation ( Documentation(revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"));
      end AlbuminTitration1;

      model AlbuminTitration2 "Figge-Fencl model (22. Dec. 2007)"
        extends Modelica.Icons.Example;

        Chemical.Solution solution(redeclare package stateOfMatter =
              Chemical.Interfaces.Incompressible)
          annotation (Placement(transformation(extent={{-104,-100},{96,100}})));

        constant Integer n=2 "Number of weak acid group in albumin molecule";
        constant Real pKAs[n]={8.5,4}
          "acid dissociation constants";
                                      //cat(1,{8.5},fill(4.0,98),fill(11.7,18),fill(12.5,24),fill(5.8,2),fill(6.0,2),{7.6,7.8,7.8,8,8},fill(10.3,50),{7.19,7.29,7.17,7.56,7.08,7.38,6.82,6.43,4.92,5.83,6.24,6.8,5.89,5.2,6.8,5.5,8,3.1})
        constant Real K[n]=fill(10.0, n) .^ (-pKAs);
        constant Real DfG[n]= Modelica.Constants.R*(298.15)*log(K);

        Chemical.Boundaries.Substance A[n](
          substanceData(each z=-1),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "deprotonated acid groups" annotation (Placement(transformation(extent={{26,-16},{6,4}})));
        Chemical.Processes.Reaction react[n](
          each nS=2,
          each nP=1) annotation (Placement(transformation(extent={{-24,-2},{-44,18}})));

        Chemical.Boundaries.Substance HA[n](
          useInlet=true,
          substanceData(DfG=DfG),
          each use_mass_start=false,
          each amountOfSubstance_start=0.00033) "protonated acid groups" annotation (Placement(transformation(extent={{-58,-2},{-78,18}})));

        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={62,-68})));
        Chemical.Boundaries.ExternalMoleFraction H(substanceData=Chemical.SubstancesOld.Proton_aqueous(), MoleFraction=10^(-7.4))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,42})));
      equation

        for i in 1:n loop
          connect(HA[i].solution, solution.solution) annotation (Line(
            points={{-62,-2},{-62,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(A[i].solution, solution.solution) annotation (Line(
            points={{22,-16},{22,-86},{56,-86},{56,-98}},
            color={127,127,0}));
          connect(react[i].substrates[1], A[i].outlet) annotation (Line(
            points={{-24,7.75},{-12,7.75},{-12,-6},{6,-6}},
            color={107,45,134},
            thickness=1));
          connect(H.outlet, react[i].substrates[2]) annotation (Line(points={{10,42},{-8,42},{-8,8.25},{-24,8.25}},
                                color={158,66,200}));
          connect(HA[i].inlet, react[i].products[1]) annotation (Line(
            points={{-58,8},{-50,8},{-50,8},{-44,8}},
            color={107,45,134},
            thickness=1));
        end for;

        connect(solution.solution, H2O.solution) annotation (Line(
          points={{56,-98},{56,-78}},
          color={127,127,0}));

        connect(H.solution, solution.solution) annotation (Line(points={{26,52},{26,78},{-94,78},{-94,-86},{56,-86},{56,-98}}, color={127,127,0}));
        annotation ( Documentation(revisions="<html>
<p><i>2014-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>",       info="<html>
<p>The titration slope der(pH)/der(SID)=185 1/(mol/L) at pH=7.4 and tAlb=0.66 mmol/l.</p>
<p>Data and model is described in</p>
<p><font style=\"color: #222222; \">Jame Figge: Role of non-volatile weak acids (albumin, phosphate and citrate). In: Stewart&apos;s Textbook of Acid-Base, 2nd Edition, John A. Kellum, Paul WG Elbers editors, &nbsp;AcidBase org, 2009, pp. 216-232.</font></p>
</html>"));
      end AlbuminTitration2;

    end AcidBase;

    package Hemoglobin "Hemoglobin blood gases binding"
      model Allosteric_Hemoglobin_MWC "Monod,Wyman,Changeux (1965)"
        extends Modelica.Icons.Example;

        constant Modelica.Units.SI.AmountOfSubstance THb=0.001
          "Total amount of hemoglobin";

        constant Modelica.Units.SI.Temperature T=298.15 "Base Temperature";
        constant Real RT=Modelica.Constants.R*T;

        constant Modelica.Units.SI.Volume OneLiter=0.001;

        constant Real L=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        constant Real c=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        constant Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
            38.7)
          "Oxygen dissociation coefficient on relaxed(R) hemoglobin subunit";

        constant Real KRx=(KR*OneLiter)
          "Mole fraction based KR";

        //Relative Gibbs formation energies of the substances in the system:
        constant Modelica.Units.SI.MolarEnergy GO2aq=-RT*log(0.0013);
        constant Modelica.Units.SI.MolarEnergy GR0=0;
        constant Modelica.Units.SI.MolarEnergy GT0=GR0 - RT*log(L);
        constant Modelica.Units.SI.MolarEnergy GR1=GR0 + GO2aq + RT*log(KRx
            /4);
        constant Modelica.Units.SI.MolarEnergy GT1=GR1 - RT*log(c*L);
        constant Modelica.Units.SI.MolarEnergy GR2=GR1 + GO2aq + RT*log(KRx
            /(3/2));
        constant Modelica.Units.SI.MolarEnergy GT2=GR2 - RT*log(c^2*L);
        constant Modelica.Units.SI.MolarEnergy GR3=GR2 + GO2aq + RT*log(KRx
            /(2/3));
        constant Modelica.Units.SI.MolarEnergy GT3=GR3 - RT*log(c^3*L);
        constant Modelica.Units.SI.MolarEnergy GR4=GR3 + GO2aq + RT*log(KRx
            *4);
        constant Modelica.Units.SI.MolarEnergy GT4=GR4 - RT*log(c^4*L);
                                        //*0.018),
        parameter Real KC = 0.0001 "Slow down factor";

        Chemical.Solution solution annotation (Placement(transformation(extent={{-68,-100},{102,106}})));

        Chemical.Boundaries.Substance oxygen_unbound(
          useInlet=true,
          substanceData(DfG=GO2aq),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-88,-70},{-68,-50}})));

        Chemical.Boundaries.Substance T0(
          useInlet=true,
          substanceData(DfG=GT0),
          use_mass_start=false,
          amountOfSubstance_start=(THb)) annotation (Placement(transformation(extent={{84,82},{64,102}})));

        Chemical.Boundaries.Substance T1(
          useInlet=true,
          substanceData(DfG=GT1),
          use_mass_start=false,
          amountOfSubstance_start=(THb*1e-4)) annotation (Placement(transformation(extent={{86,42},{66,62}})));

        Chemical.Boundaries.Substance T2(
          useInlet=true,
          substanceData(DfG=GT2),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-8) annotation (Placement(transformation(extent={{86,6},{66,26}})));

        Chemical.Boundaries.Substance R1(
          useInlet=true,
          substanceData(DfG=GR1),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-8) annotation (Placement(transformation(extent={{-26,42},{-46,62}})));

        Chemical.Boundaries.Substance R2(
          useInlet=true,
          substanceData(DfG=GR2),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-10) annotation (Placement(transformation(extent={{-26,6},{-46,26}})));

        Chemical.Boundaries.Substance T3(
          useInlet=true,
          substanceData(DfG=GT3),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-12) annotation (Placement(transformation(extent={{88,-56},{68,-36}})));

        Chemical.Boundaries.Substance R3(
          useInlet=true,
          substanceData(DfG=GR3),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-12) annotation (Placement(transformation(extent={{-26,-56},{-46,-36}})));

        Chemical.Boundaries.Substance T4(
          useInlet=false,
          substanceData(DfG=GT4),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-17) annotation (Placement(transformation(extent={{88,-92},{68,-72}})));

        Chemical.Boundaries.Substance R4(
          useInlet=true,
          substanceData(DfG=GR4),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-14) annotation (Placement(transformation(extent={{-26,-92},{-46,-72}})));

        Chemical.Boundaries.Substance R0(
          useInlet=true,
          substanceData(DfG=GR0),
          use_mass_start=false,
          amountOfSubstance_start=THb*1e-7) annotation (Placement(transformation(extent={{-26,82},{-46,102}})));

        Chemical.Processes.Reaction quaternaryForm(
          nP=1,
          nS=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,92})));

        Chemical.Processes.Reaction oxyR1(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-54,66})));

        Chemical.Processes.Reaction oxyT1(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={70,72})));

        Chemical.Processes.Reaction oxyR2(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-54,28})));

        Chemical.Processes.Reaction oxyR3(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-54,-26})));

        Chemical.Processes.Reaction oxyR4(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-52,-70})));

        Chemical.Processes.Reaction oxyT2(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={72,34})));

        Chemical.Processes.Reaction oxyT3(
          nS=1,
          nP=2)  annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={74,-28})));

        Chemical.Processes.Reaction oxyT4(
          nP=2,
          nS=1)  annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={74,-64})));

        Chemical.Processes.Reaction quaternaryForm1(
          nS=1,
          nP=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,52})));

        Chemical.Processes.Reaction quaternaryForm2(
          nP=1,
          nS=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,16})));

        Chemical.Processes.Reaction quaternaryForm3(
          nP=1,
          nS=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,-46})));

        Chemical.Processes.Reaction quaternaryForm4(
          nP=1,
          nS=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={22,-82})));

        Modelica.Blocks.Sources.ContinuousClock clock(offset=10)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-82,78})));
        Chemical.Boundaries.ExternalIdealGas O2_in_air(
          usePartialPressureInput=true,
          substanceData=Chemical.SubstancesOld.Oxygen_gas(),
          PartialPressure(displayUnit="kPa") = 3733)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-82,42})));

        Chemical.Processes.GasSolubility gasSolubility annotation (Placement(transformation(extent={{-92,0},{-72,20}})));

        Real sO2;
        Chemical.Boundaries.Substance substance(substanceData=SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{66,-20},{86,0}})));
        Chemical.Topology.SplitterT2 splitterT2 annotation (Placement(transformation(extent={{58,42},{38,62}})));
        TopologyToProcess.JunctionT2 junctionT2 annotation (Placement(transformation(extent={{-18,102},{2,82}})));
        TopologyToProcess.JunctionT2 junctionT2_1 annotation (Placement(transformation(extent={{-18,62},{2,42}})));
        TopologyToProcess.JunctionT2 junctionT2_2 annotation (Placement(transformation(extent={{-18,26},{2,6}})));
        Chemical.Topology.SplitterT2 splitterT1 annotation (Placement(transformation(extent={{58,6},{38,26}})));
        Chemical.Topology.SplitterT2 splitterT3 annotation (Placement(transformation(extent={{58,-56},{38,-36}})));
        TopologyToProcess.JunctionT2 junctionT2_3 annotation (Placement(transformation(extent={{-18,-36},{2,-56}})));
        Chemical.Topology.SplitterT2 splitterT4 annotation (Placement(transformation(extent={{60,-92},{40,-72}})));
        TopologyToProcess.JunctionN junctionN(N=9)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-82,-16})));
        inner DropOfCommons dropOfCommons(L=10) annotation (Placement(transformation(extent={{-112,-92},{-92,-72}})));
      equation
        sO2 = (R1.substance.x + 2*R2.substance.x + 3*R3.substance.x + 4*R4.substance.x + T1.substance.x + 2*T2.substance.x + 3*T3.substance.x + 4*T4.substance.x) /
         (4*(R0.substance.x + R1.substance.x + R2.substance.x + R3.substance.x + R4.substance.x + T0.substance.x + T1.substance.x + T2.substance.x + T3.substance.x + T4.substance.x));

        connect(oxygen_unbound.solution, solution.solution) annotation (Line(
            points={{-84,-70},{-84,-112},{102,-112},{102,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}},
            color={127,127,0}));
        connect(T0.solution, solution.solution) annotation (Line(
            points={{80,82},{106,82},{106,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}},
            color={127,127,0}));
        connect(T1.solution, solution.solution) annotation (Line(points={{82,42},{106,42},{106,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}},
                            color={127,127,0}));
        connect(T2.solution, solution.solution) annotation (Line(points={{82,6},{106,6},{106,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}},
                                  color={127,127,0}));
        connect(substance.solution, solution.solution) annotation (Line(points={{70,-20},{106,-20},{106,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}},
                                                             color={127,127,0}));
        connect(clock.y, O2_in_air.partialPressure)
          annotation (Line(points={{-82,67},{-82,52}},          color={0,0,127}));
        connect(O2_in_air.solution, solution.solution) annotation (Line(points={{-92,48},{-96,48},{-96,-56},{-66,-56},{-66,-104},{-30,-104},{-30,-108},{68,-108},{68,-97.94}},
                                                                                                                                 color={127,127,0}));
        connect(O2_in_air.outlet, gasSolubility.inlet) annotation (Line(
            points={{-82,32},{-82,22},{-82,10},{-92,10}},
            color={158,66,200},
            thickness=0.5));
        connect(R0.solution, solution.solution)
          annotation (Line(points={{-30,82},{-30,110},{-66,110},{-66,-104},{-30,-104},{-30,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(R1.solution, solution.solution)
          annotation (Line(points={{-30,42},{-66,42},{-66,-104},{-30,-104},{-30,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(R2.solution, solution.solution) annotation (Line(points={{-30,6},{-66,6},{-66,-104},{-30,-104},{-30,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(junctionT2_1.outlet, R1.inlet) annotation (Line(
            points={{-18,52},{-26,52}},
            color={158,66,200},
            thickness=0.5));
        connect(R3.solution, solution.solution)
          annotation (Line(points={{-30,-56},{-66,-56},{-66,-104},{-30,-104},{-30,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(R4.solution, solution.solution) annotation (Line(points={{-30,-92},{-30,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(T3.solution, solution.solution)
          annotation (Line(points={{84,-56},{106,-56},{106,-104},{84,-104},{84,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(T4.solution, solution.solution) annotation (Line(points={{84,-92},{84,-108},{68,-108},{68,-97.94}}, color={127,127,0}));
        connect(junctionT2.inletB, quaternaryForm.products[1]) annotation (Line(
            points={{2,92},{10,92}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm.substrates[1], T0.outlet) annotation (Line(
            points={{30,92},{64,92}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT2.inletProcess, T1.outlet) annotation (Line(
            points={{58,52},{66,52}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm1.substrates[1], splitterT2.outletSubstance1)
          annotation (Line(
            points={{30,52},{38,52}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT2_1.inletB, quaternaryForm1.products[1]) annotation (Line(
            points={{2,52},{10,52}},
            color={158,66,200},
            thickness=0.5));
        connect(R0.inlet, junctionT2.outlet) annotation (Line(
            points={{-26,92},{-18,92}},
            color={158,66,200},
            thickness=0.5));
        connect(R2.inlet, junctionT2_2.outlet) annotation (Line(
            points={{-26,16},{-18,16}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT2_2.inletB, quaternaryForm2.products[1]) annotation (Line(
            points={{2,16},{10,16}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm2.substrates[1], splitterT1.outletSubstance1)
          annotation (Line(
            points={{30,16},{38,16}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.inletProcess, T2.outlet) annotation (Line(
            points={{58,16},{66,16}},
            color={158,66,200},
            thickness=0.5));
        connect(R3.inlet, junctionT2_3.outlet) annotation (Line(
            points={{-26,-46},{-18,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT2_3.inletB, quaternaryForm3.products[1]) annotation (Line(
            points={{2,-46},{10,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm3.substrates[1], splitterT3.outletSubstance1)
          annotation (Line(
            points={{30,-46},{38,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT3.inletProcess, T3.outlet) annotation (Line(
            points={{58,-46},{68,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm4.products[1], R4.inlet) annotation (Line(
            points={{12,-82},{-26,-82}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm4.substrates[1], splitterT4.outletSubstance1)
          annotation (Line(
            points={{32,-82},{40,-82}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT4.inletProcess, T4.outlet) annotation (Line(
            points={{60,-82},{68,-82}},
            color={158,66,200},
            thickness=0.5));
        connect(R4.outlet, oxyR4.substrates[1]) annotation (Line(
            points={{-46,-82},{-46,-80},{-52,-80}},
            color={158,66,200},
            thickness=0.5));
        connect(R3.outlet, oxyR3.substrates[1]) annotation (Line(
            points={{-46,-46},{-54,-46},{-54,-36}},
            color={158,66,200},
            thickness=0.5));
        connect(R2.outlet, oxyR2.substrates[1]) annotation (Line(
            points={{-46,16},{-46,18},{-54,18}},
            color={158,66,200},
            thickness=0.5));
        connect(R1.outlet, oxyR1.substrates[1]) annotation (Line(
            points={{-46,52},{-54,52},{-54,56}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR1.products[1], junctionT2.inletA) annotation (Line(
            points={{-54.25,76},{-54.25,78},{-8,78},{-8,82}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR2.products[1], junctionT2_1.inletA) annotation (Line(
            points={{-54.25,38},{-8,38},{-8,42}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR3.products[1], junctionT2_2.inletA) annotation (Line(
            points={{-54.25,-16},{-8,-16},{-8,6}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR4.products[1], junctionT2_3.inletA) annotation (Line(
            points={{-52.25,-60},{-8,-60},{-8,-56}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT4.outletSubstance, oxyT4.substrates[1])
          annotation (Line(
            points={{50,-72},{50,-64},{64,-64}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT3.outletSubstance, oxyT3.substrates[1])
          annotation (Line(
            points={{48,-36},{48,-28},{64,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletSubstance, oxyT2.substrates[1]) annotation (Line(
            points={{48,26},{48,34},{62,34}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT2.outletSubstance, oxyT1.substrates[1]) annotation (Line(
            points={{48,62},{48,72},{60,72}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT4.products[1], T3.inlet) annotation (Line(
            points={{84,-64.25},{90,-64.25},{90,-46},{88,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT3.products[1], T2.inlet) annotation (Line(
            points={{84,-28.25},{90,-28.25},{90,16},{86,16}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT2.products[1], T1.inlet) annotation (Line(
            points={{82,33.75},{90,33.75},{90,52},{86,52}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT1.products[1], T0.inlet) annotation (Line(
            points={{80,71.75},{90,71.75},{90,92},{84,92}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionN.outlet, oxygen_unbound.inlet) annotation (Line(
            points={{-92,-16},{-100,-16},{-100,-60},{-88,-60}},
            color={158,66,200},
            thickness=0.5));
        connect(O2_in_air.solution, oxygen_unbound.solution) annotation (Line(points={{-92,48},{-118,48},{-118,-96},{-84,-96},{-84,-70}},
                                                                                                                              color={127,127,0}));
        connect(gasSolubility.outlet, junctionN.inlets[1]) annotation (Line(
            points={{-82,0},{-72,0},{-72,-15.1111}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR4.products[2], junctionN.inlets[2])
          annotation (Line(
            points={{-51.75,-60},{-64,-60},{-64,-15.3333},{-72,-15.3333}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR3.products[2], junctionN.inlets[3])
          annotation (Line(
            points={{-53.75,-16},{-60,-16},{-60,-15.5556},{-72,-15.5556}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyR2.products[2], junctionN.inlets[4])
          annotation (Line(
            points={{-53.75,38},{-64,38},{-64,-15.7778},{-72,-15.7778}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionN.inlets[5], oxyR1.products[2])
          annotation (Line(
            points={{-72,-16},{-64,-16},{-64,78},{-53.75,78},{-53.75,76}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT4.products[2], junctionN.inlets[6])
          annotation (Line(
            points={{84,-63.75},{100,-63.75},{100,-110},{-72,-110},{-72,-16.2222}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT3.products[2], junctionN.inlets[7])
          annotation (Line(
            points={{84,-27.75},{100,-27.75},{100,-110},{-72,-110},{-72,-16.4444}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT2.products[2], junctionN.inlets[8])
          annotation (Line(
            points={{82,34.25},{100,34.25},{100,-110},{-72,-110},{-72,-16.6667}},
            color={158,66,200},
            thickness=0.5));
        connect(oxyT1.products[2], junctionN.inlets[9])
          annotation (Line(
            points={{80,72.25},{100,72.25},{100,-110},{-72,-110},{-72,-16.8889}},
            color={158,66,200},
            thickness=0.5));
        annotation ( Documentation(info="<html>
<p>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &quot;On the nature of allosteric transitions: a plausible model.&quot; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 %).</p>
</html>",   revisions="<html>
<p><i>2013-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),   experiment(StopTime=15000, __Dymola_Algorithm="Dassl"),
          __Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin_MWC;

      model Allosteric_Hemoglobin2_MWC
        "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"
        extends Modelica.Icons.Example;

        constant Modelica.Units.SI.AmountOfSubstance THb=0.001
          "Total amount of hemoglobin";

        constant Modelica.Units.SI.Temperature T=298.15 "Base Temperature";
        constant Real RT=Modelica.Constants.R*T;

        // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        //   "Amount of solution used for molarity to mole fraction conversion";
        constant Modelica.Units.SI.Volume OneLiter=0.001;

        parameter Real L=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        parameter Real c=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
            38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                    //*7.875647668393782383419689119171e-5
        //10.500001495896 7.8756465463794e-05
        parameter Modelica.Units.SI.Concentration KT=KR/c
          "oxygen dissociation on tensed(T) hemoglobin subunit";

        parameter Modelica.Units.SI.MoleFraction KRx=KR*OneLiter;
        //AmountOfSolutionIn1L;
        parameter Modelica.Units.SI.MoleFraction KTx=KT*OneLiter;
        //AmountOfSolutionIn1L;
        parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
        parameter Modelica.Units.SI.ChemicalPotential DfG_uR=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_uRO2=DfG_uR +
            DfG_O2 + RT*log(KRx);
        parameter Modelica.Units.SI.ChemicalPotential DfG_uT=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_uTO2=DfG_uT +
            DfG_O2 + RT*log(KTx);
        parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
            log(L);

        parameter Real KC = 1e-6 "Slow down factor";
                                 //0.000001
        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        Chemical.Processes.Reaction quaternaryForm(
          nS=1,
          nP=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={22,-48})));
        Chemical.Processes.SpeciationIn R0_in_R(NumberOfSubunits=4) annotation (Placement(transformation(extent={{-46,-48},{-26,-28}})));
         // AmountOfSubstance_start=4e-11)
        Chemical.Processes.SpeciationOut T0_in_T(NumberOfSubunits=4) annotation (Placement(transformation(extent={{74,-48},{54,-28}})));
         // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Chemical.Boundaries.Substance OxyRHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KRx) + DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=2e-12) "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-76,-20},{-96,0}})));

        Chemical.Processes.Reaction oxygenation_R[4](
          each nS=2, each nP=1)             annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-52,-10})));
        Chemical.Boundaries.Substance DeoxyRHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1.471e-10) "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-8,-20},{-28,0}})));

        Chemical.Boundaries.Substance OxyTHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KTx) + DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=5e-8) "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{28,-16},{8,4}})));

        Chemical.Processes.Reaction oxygenation_T[4](
          each nS=2, each nP=1)
                      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={42,-8})));
        Chemical.Boundaries.Substance DeoxyTHm[4](
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1e-3) "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{106,-18},{86,2}})));

        Chemical.Boundaries.Substance oxygen_unbound(
          useInlet=true,
          substanceData(DfG=DfG_O2),
          use_mass_start=false,
          amountOfSubstance_start=2e-9) annotation (Placement(transformation(extent={{2,34},{22,54}})));
        Modelica.Blocks.Sources.ContinuousClock clock(offset=1) annotation (
           Placement(transformation(extent={{-82,70},{-62,90}})));
        Chemical.Boundaries.ExternalIdealGas oxygen_in_air(usePartialPressureInput=true)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,80})));
        Chemical.Processes.GasSolubility partialPressure1 annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,62})));

        Real sO2 "Hemoglobin oxygen saturation";
        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{32,-92},{52,-72}})));
        Chemical.Topology.SplitterT2 sT2[4] annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={66,-8})));
        Chemical.Topology.SplitterN uO2(N=8) annotation (Placement(transformation(extent={{42,34},{62,54}})));
        inner DropOfCommons dropOfCommons(L=1e-5) annotation (Placement(transformation(extent={{56,68},{82,92}})));
      equation
        sO2 = (sum(OxyRHm.x) + sum(OxyTHm.x)) /
        (sum(DeoxyRHm.x) + sum(DeoxyTHm.x) + sum(OxyRHm.x) + sum(OxyTHm.x));

        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-61,80},{-46,80}},
            color={0,0,127}));

        for i in 1:4 loop

          connect(oxygenation_T[i].substrates[2], uO2.outlet[i]) annotation (Line(points={{52,-8.25},{62,-8.25},{62,44}}, color={107,45,134}));
          connect(oxygenation_R[i].substrates[2], uO2.outlet[i + 4])
            annotation (Line(points={{-42,-10.25},{-38,-10.25},{-38,20},{62,20},{62,44}}, color={107,45,134}));

        connect(R0_in_R.subunitSolution, DeoxyRHm[i].solution) annotation (Line(
            points={{-32,-32},{-32,-22},{-12,-22},{-12,-20}},
            color={127,127,0}));
        connect(R0_in_R.subunitSolution, OxyRHm[i].solution) annotation (Line(
            points={{-32,-32},{-32,-22},{-80,-22},{-80,-20}},
            color={127,127,0}));
        connect(OxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{24,-16},{24,-22},{60,-22},{60,-32}},
            color={127,127,0}));
        connect(DeoxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{102,-18},{102,-22},{60,-22},{60,-32}},
            color={127,127,0}));
        end for;

        connect(R0_in_R.solution, solution.solution) annotation (Line(
            points={{-42,-48},{0,-48},{0,-98},{60,-98}},
            color={127,127,0}));
        connect(T0_in_T.solution, solution.solution) annotation (Line(
            points={{70,-48},{104,-48},{104,-104},{60,-104},{60,-98}},
            color={127,127,0}));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(points={{6,34},{6,-100},{104,-100},{104,-104},{60,-104},{60,-98}},
                                                 color={127,127,0}));
        connect(solution.solution, H2O.solution) annotation (Line(
            points={{60,-98},{36,-98},{36,-92}},
            color={127,127,0}));

        connect(T0_in_T.outletSubstance, quaternaryForm.substrates[1])
          annotation (Line(
            points={{54,-48},{32,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm.products[1], R0_in_R.inlet) annotation (Line(
            points={{12,-48},{-26,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_in_air.outlet, partialPressure1.inlet) annotation (Line(
            points={{-26,80},{-22,80},{-22,62},{-24,62}},
            color={158,66,200},
            thickness=0.5));
        connect(partialPressure1.outlet, oxygen_unbound.inlet)
          annotation (Line(
            points={{-14,52},{-14,42},{2,42},{2,44}},
            color={158,66,200},
            thickness=0.5));
        connect(R0_in_R.subunits, DeoxyRHm.inlet)
          annotation (Line(
            points={{-39,-27.8},{-4,-27.8},{-4,-10},{-8,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_R.substrates[1], DeoxyRHm.outlet) annotation (Line(
            points={{-42,-9.75},{-36,-9.75},{-36,-10},{-28,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_R.products[1], OxyRHm.inlet) annotation (Line(
            points={{-62,-10},{-76,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_T.products[1], OxyTHm.inlet) annotation (Line(
            points={{32,-8},{32,-6},{28,-6}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_in_air.solution, solution.solution) annotation (Line(points={{-42,70},{-42,46},{-98,46},{-98,-156},{60,-156},{60,-98}}, color={127,127,0}));
        connect(DeoxyTHm.outlet, sT2.inletProcess) annotation (Line(
            points={{86,-8},{76,-8}},
            color={158,66,200},
            thickness=0.5));
        connect(sT2.outletSubstance, T0_in_T.subunits)
          annotation (Line(
            points={{66,-18},{66,-22.9},{67,-22.9},{67,-27.8}},
            color={158,66,200},
            thickness=0.5));
        connect(sT2.outletSubstance1, oxygenation_T.substrates[1])
          annotation (Line(
            points={{56,-8},{54,-8},{54,-7.75},{52,-7.75}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_unbound.outlet, uO2.inlet) annotation (Line(
            points={{22,44},{42,44}},
            color={158,66,200},
            thickness=0.5));
        annotation (          experiment(StopTime=15000, __Dymola_Algorithm="Dassl"),
          Documentation(revisions="<html>
<p><i>2013-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p><br>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &quot;On the nature of allosteric transitions: a plausible model.&quot; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 %).</p>
</html>"),__Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin2_MWC;

      model Allosteric_Hemoglobin2_MWC_
        "Monod,Wyman,Changeux (1965) - The same allosteric hemoglobin model as Allosteric_Hemoglobin_MWC implemented by Speciation blocks"
        extends Modelica.Icons.Example;

        constant Modelica.Units.SI.AmountOfSubstance THb=0.001
          "Total amount of hemoglobin";

        constant Modelica.Units.SI.Temperature T=298.15 "Base Temperature";
        constant Real RT=Modelica.Constants.R*T;

        // constant Modelica.SIunits.AmountOfSubstance AmountOfSolutionIn1L = 38.7
        //   "Amount of solution used for molarity to mole fraction conversion";
        constant Modelica.Units.SI.Volume OneLiter=0.001;

        parameter Real L=7.0529*10^6
          "=[T0]/[R0] .. dissociation constant of relaxed <-> tensed change of deoxyhemoglobin tetramer";
        parameter Real c=0.00431555
          "=KR/KT .. ration between oxygen affinities of relaxed vs. tensed subunit";
        parameter Modelica.Units.SI.Concentration KR=0.000671946*(55.508/
            38.7) "oxygen dissociation on relaxed(R) hemoglobin subunit";
                                                                    //*7.875647668393782383419689119171e-5
        //10.500001495896 7.8756465463794e-05
        parameter Modelica.Units.SI.Concentration KT=KR/c
          "oxygen dissociation on tensed(T) hemoglobin subunit";

        parameter Modelica.Units.SI.MoleFraction KRx=KR*OneLiter;
        //AmountOfSolutionIn1L;
        parameter Modelica.Units.SI.MoleFraction KTx=KT*OneLiter;
        //AmountOfSolutionIn1L;
        parameter Modelica.Units.SI.ChemicalPotential DfG_O2=-RT*log(0.0013);
        parameter Modelica.Units.SI.ChemicalPotential DfG_uR=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_uRO2=DfG_uR +
            DfG_O2 + RT*log(KRx);
        parameter Modelica.Units.SI.ChemicalPotential DfG_uT=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_uTO2=DfG_uT +
            DfG_O2 + RT*log(KTx);
        parameter Modelica.Units.SI.ChemicalPotential DfG_tT=0;
        parameter Modelica.Units.SI.ChemicalPotential DfG_tR=DfG_tT + RT*
            log(L);

        parameter Real KC = 1e-6 "Slow down factor";
                                 //0.000001
        Chemical.Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

        Chemical.Processes.Reaction quaternaryForm(
          nS=1,
          nP=1)  annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={22,-48})));
        Chemical.Processes.SpeciationIn R0_in_R(NumberOfSubunits=4) annotation (Placement(transformation(extent={{-46,-48},{-26,-28}})));
         // AmountOfSubstance_start=4e-11)
        Chemical.Processes.SpeciationOut T0_in_T(NumberOfSubunits=4) annotation (Placement(transformation(extent={{74,-48},{54,-28}})));
         // AmountOfSubstance_start=totalAmountOfHemoglobin)
        Chemical.Boundaries.Substance OxyRHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KRx) + DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=2e-12) "Oxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-76,-20},{-96,0}})));

        Chemical.Processes.Reaction oxygenation_R[4](
          each nS=2, each nP=1)             annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-52,-10})));
        Chemical.Boundaries.Substance DeoxyRHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tR/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1.471e-10) "Deoxygenated subunit in R structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{-8,-20},{-28,0}})));

        Chemical.Boundaries.Substance OxyTHm[4](
          each useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_O2 + RT*log(KTx) + DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=5e-8) "Oxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{28,-16},{8,4}})));

        Chemical.Processes.Reaction oxygenation_T[4](
          each nS=2, each nP=1)
                      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=180,
              origin={42,-8})));
        Chemical.Boundaries.Substance DeoxyTHm[4](
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          each substanceData(DfG=DfG_tT/4),
          each use_mass_start=false,
          each amountOfSubstance_start=1e-3) "Deoxygenated subunit in T structure of hemoglobin tetramer"
          annotation (Placement(transformation(extent={{106,-18},{86,2}})));

        Chemical.Boundaries.Substance oxygen_unbound(
          useInlet=true,
          substanceData(DfG=DfG_O2),
          use_mass_start=false,
          amountOfSubstance_start=2e-9) annotation (Placement(transformation(extent={{2,34},{22,54}})));
        Modelica.Blocks.Sources.ContinuousClock clock(offset=1) annotation (
           Placement(transformation(extent={{-82,70},{-62,90}})));
        Chemical.Boundaries.ExternalIdealGas oxygen_in_air(usePartialPressureInput=true)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,80})));
        Chemical.Processes.GasSolubility partialPressure1 annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={-14,62})));

        Real sO2 "Hemoglobin oxygen saturation";
        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{32,-92},{52,-72}})));
        Chemical.Topology.SplitterT2 sT2[4] annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={66,-8})));
        Chemical.Topology.SplitterN uO2(N=8) annotation (Placement(transformation(extent={{42,34},{62,54}})));
        inner DropOfCommons dropOfCommons(L=1e-5) annotation (Placement(transformation(extent={{56,68},{82,92}})));
      equation
        sO2 = (sum(OxyRHm.x) + sum(OxyTHm.x)) /
        (sum(DeoxyRHm.x) + sum(DeoxyTHm.x) + sum(OxyRHm.x) + sum(OxyTHm.x));

        connect(clock.y, oxygen_in_air.partialPressure) annotation (Line(
            points={{-61,80},{-46,80}},
            color={0,0,127}));

        for i in 1:4 loop

          connect(oxygenation_T[i].substrates[2], uO2.outlet[i]) annotation (Line(points={{52,-8.25},{62,-8.25},{62,44}}, color={107,45,134}));
          connect(oxygenation_R[i].substrates[2], uO2.outlet[i + 4])
            annotation (Line(points={{-42,-10.25},{-38,-10.25},{-38,20},{62,20},{62,44}}, color={107,45,134}));

        connect(R0_in_R.subunitSolution, DeoxyRHm[i].solution) annotation (Line(
            points={{-32,-32},{-32,-104},{-12,-104},{-12,-20}},
            color={127,127,0}));
        connect(R0_in_R.subunitSolution, OxyRHm[i].solution) annotation (Line(
            points={{-32,-32},{-32,-22},{-80,-22},{-80,-20}},
            color={127,127,0}));
        connect(OxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{24,-16},{24,-22},{60,-22},{60,-32}},
            color={127,127,0}));
        connect(DeoxyTHm[i].solution, T0_in_T.subunitSolution) annotation (Line(
            points={{102,-18},{102,-22},{60,-22},{60,-32}},
            color={127,127,0}));
        end for;

        connect(R0_in_R.solution, solution.solution) annotation (Line(
            points={{-42,-48},{0,-48},{0,-98},{60,-98}},
            color={127,127,0}));
        connect(T0_in_T.solution, solution.solution) annotation (Line(
            points={{70,-48},{104,-48},{104,-104},{60,-104},{60,-98}},
            color={127,127,0}));
        connect(oxygen_unbound.solution, solution.solution) annotation (Line(points={{6,34},{6,-100},{104,-100},{104,-104},{60,-104},{60,-98}},
                                                 color={127,127,0}));
        connect(solution.solution, H2O.solution) annotation (Line(
            points={{60,-98},{36,-98},{36,-92}},
            color={127,127,0}));

        connect(T0_in_T.outletSubstance, quaternaryForm.substrates[1])
          annotation (Line(
            points={{54,-48},{32,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(quaternaryForm.products[1], R0_in_R.inlet) annotation (Line(
            points={{12,-48},{-26,-48}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_in_air.outlet, partialPressure1.inlet) annotation (Line(
            points={{-26,80},{-22,80},{-22,72},{-14,72}},
            color={158,66,200},
            thickness=0.5));
        connect(partialPressure1.outlet, oxygen_unbound.inlet)
          annotation (Line(
            points={{-14,52},{-14,42},{2,42},{2,44}},
            color={158,66,200},
            thickness=0.5));
        connect(R0_in_R.subunits, DeoxyRHm.inlet)
          annotation (Line(
            points={{-39,-27.8},{-8,-27.8},{-8,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_R.substrates[1], DeoxyRHm.outlet) annotation (Line(
            points={{-42,-9.75},{-42,-10},{-28,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_R.products[1], OxyRHm.inlet) annotation (Line(
            points={{-62,-10},{-76,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygenation_T.products[1], OxyTHm.inlet) annotation (Line(
            points={{32,-8},{32,-6},{28,-6}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_in_air.solution, solution.solution) annotation (Line(points={{-42,70},{-42,46},{-98,46},{-98,-156},{60,-156},{60,-98}}, color={127,127,0}));
        connect(DeoxyTHm.outlet, sT2.inletProcess) annotation (Line(
            points={{86,-8},{76,-8}},
            color={158,66,200},
            thickness=0.5));
        connect(sT2.outletSubstance, T0_in_T.subunits)
          annotation (Line(
            points={{66,-18},{66,-22.9},{67,-22.9},{67,-27.8}},
            color={158,66,200},
            thickness=0.5));
        connect(sT2.outletSubstance1, oxygenation_T.substrates[1])
          annotation (Line(
            points={{56,-8},{54,-8},{54,-7.75},{52,-7.75}},
            color={158,66,200},
            thickness=0.5));
        connect(oxygen_unbound.outlet, uO2.inlet) annotation (Line(
            points={{22,44},{42,44}},
            color={158,66,200},
            thickness=0.5));
        annotation (          experiment(StopTime=15000, __Dymola_Algorithm="Dassl"),
          Documentation(revisions="<html>
<p><i>2013-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>", info="<html>
<p><br>To understand the model is necessary to study the principles of MWC allosteric transitions first published by </p>
<p>[1] Monod,Wyman,Changeux (1965). &quot;On the nature of allosteric transitions: a plausible model.&quot; Journal of molecular biology 12(1): 88-118.</p>
<p><br>In short it is about binding oxygen to hemoglobin.</p>
<p>Oxgen are driven by its partial pressure using clock source - from very little pressure to pressure of 10kPa.</p>
<p>(Partial pressure of oxygen in air is the air pressure multiplied by the fraction of the oxygen in air.)</p>
<p>Hemoglobin was observed (by Perutz) in two structuraly different forms R and T.</p>
<p>These forms are represented by blocks T0..T4 and R0..R4, where the suffexed index means the number of oxygen bounded to the form.</p>
<p><br>In equilibrated model can be four chemical reactions removed and the results will be the same, but dynamics will change a lot. ;)</p>
<p>If you remove the quaternaryForm1,quaternaryForm2,quaternaryForm3,quaternaryForm4 then the model in equilibrium will be exactly the same as in MWC article.</p>
<p><br>Parameters was fitted to data of Severinghaus article from 1979. (For example at pO2=26mmHg is oxygen saturation sO2 = 48.27 %).</p>
</html>"),__Dymola_experimentSetupOutput);
      end Allosteric_Hemoglobin2_MWC_;
    end Hemoglobin;

    model H2O_ElectrochemicalCell
     extends Modelica.Icons.Example;

      Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
      Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

      Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-96},{56,-30}})));

      Chemical.Boundaries.ExternalIdealGas H2(substanceData=Chemical.SubstancesOld.Hydrogen_gas(), PartialPressure=100000)
        annotation (Placement(transformation(extent={{44,32},{24,52}})));
      Chemical.Boundaries.Substance H(
        substanceData=Chemical.SubstancesOld.Proton_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-4,-58},{16,-38}})));
      Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
        annotation (Placement(transformation(extent={{-6,64},{14,84}})));
      Chemical.Processes.Reaction electrodeReaction(
        s={2,2},
        nS=2,
        nP=1) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={52,6})));
      Chemical.Processes.Reaction electrodeReaction1(
        s={2,8,8},
        p={4},
        nP=1,
        nS=3) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={-40,0})));

      Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                 //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                  //(substanceData=Chemical.Examples.Substances.Electrone_solid())
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{62,54},{82,74}})));
      Chemical.Boundaries.Substance H2O(
        useInlet=true,
        substanceData=Chemical.SubstancesOld.Water_liquid(),
        mass_start=1) annotation (Placement(transformation(extent={{-6,-82},{14,-62}})));
      Chemical.Boundaries.ExternalIdealGas O2_(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure=100000)
        annotation (Placement(transformation(extent={{-4,36},{-24,56}})));
      Solution gases(BasePressure(displayUnit="kPa") = 100000) annotation (Placement(transformation(extent={{-42,22},{60,58}})));
      Chemical.Topology.SplitterT1 splitterT1 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={20,6})));
      Chemical.Boundaries.Substance O2_aq(
        useInlet=true,
        substanceData=Chemical.SubstancesOld.Oxygen_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 0.0001) annotation (Placement(transformation(extent={{-28,-62},{-8,-42}})));
      Chemical.Boundaries.Substance H2_aq(
        useInlet=true,
        substanceData=Chemical.SubstancesOld.Hydrogen_aqueous(),
        use_mass_start=false,
        amountOfSubstance_start(displayUnit="mmol") = 1e-05) annotation (Placement(transformation(extent={{24,-76},{44,-56}})));
      Chemical.Processes.GasSolubility gasSolubility annotation (Placement(transformation(extent={{-36,-22},{-16,-2}})));
      Chemical.Processes.GasVolatility gasVolatility annotation (Placement(transformation(extent={{40,-46},{60,-26}})));
    equation
      connect(H.solution, solution1.solution) annotation (Line(points={{0,-58},{0,-26},{-42,-26},{-42,-102},{-2,-102},{-2,-100},{38.8,-100},{38.8,-95.34}},
                                         color={127,127,0}));
    connect(electrone.solution, cathode.solution) annotation (Line(
        points={{-74,32},{-74,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
        color={127,127,0}));
    connect(electrone1.solution, anode.solution) annotation (Line(
        points={{84,-26},{84,-49},{89.2,-49}},
        color={127,127,0}));
      connect(voltageSensor.p, electrone.pin) annotation (Line(
          points={{-6,74},{-96,74},{-96,51.8},{-68,51.8}},
          color={0,0,255}));
      connect(voltageSensor.n, electrone1.pin) annotation (Line(
          points={{14,74},{92,74},{92,-6.2},{78,-6.2}},
          color={0,0,255}));
      connect(electrone1.pin, ground.p) annotation (Line(
          points={{78,-6.2},{92,-6.2},{92,74},{72,74}},
          color={0,0,255}));
      connect(H2O.solution, solution1.solution) annotation (Line(points={{-2,-82},{-2,-100},{38.8,-100},{38.8,-95.34}},
                                          color={127,127,0}));
      connect(gases.solution, H2.solution) annotation (Line(points={{39.6,22.36},{40,22.36},{40,32}}, color={127,127,0}));
      connect(O2_.solution, gases.solution) annotation (Line(points={{-8,36},{-8,18},{39.6,18},{39.6,22.36}}, color={127,127,0}));
      connect(electrodeReaction1.products[1], H2O.inlet) annotation (Line(
          points={{-40,-10},{-40,-72},{-6,-72}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterT1.outletB, electrodeReaction.substrates[1])
        annotation (Line(
          points={{30,6},{38,6},{38,-16},{51.75,-16},{51.75,-4}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone1.outlet, electrodeReaction.substrates[2])
        annotation (Line(
          points={{68,-16},{52.25,-16},{52.25,-4}},
          color={158,66,200},
          thickness=0.5));
      connect(H.outlet, splitterT1.inlet) annotation (Line(
          points={{16,-48},{16,-12},{20,-12},{20,-4}},
          color={158,66,200},
          thickness=0.5));
      connect(O2_aq.solution, solution1.solution) annotation (Line(points={{-24,-62},{-24,-95.34},{38.8,-95.34}}, color={127,127,0}));
      connect(H2_aq.solution, solution1.solution) annotation (Line(points={{28,-76},{28,-100},{38.8,-100},{38.8,-95.34}}, color={127,127,0}));
      connect(O2_aq.outlet, electrodeReaction1.substrates[1])
        annotation (Line(
          points={{-8,-52},{-16,-52},{-16,16},{-40.3333,16},{-40.3333,10}},
          color={158,66,200},
          thickness=0.5));
      connect(electrone.outlet, electrodeReaction1.substrates[2])
        annotation (Line(
          points={{-58,42},{-40,42},{-40,14},{-40,14},{-40,10}},
          color={158,66,200},
          thickness=0.5));
      connect(splitterT1.outletA, electrodeReaction1.substrates[3])
        annotation (Line(
          points={{10,6},{-26,6},{-26,18},{-39.6667,18},{-39.6667,10}},
          color={158,66,200},
          thickness=0.5));
      connect(gasSolubility.outlet, O2_aq.inlet) annotation (Line(
          points={{-16,-12},{-16,-52},{-28,-52}},
          color={200,66,175},
          thickness=0.5));
      connect(gasSolubility.inlet, O2_.outlet) annotation (Line(
          points={{-36,-12},{-36,46},{-24,46}},
          color={158,66,200},
          thickness=0.5));
      connect(electrodeReaction.products[1], H2_aq.inlet)
        annotation (Line(
          points={{52,16},{34,16},{34,-48},{20,-48},{20,-66},{24,-66}},
          color={200,66,175},
          thickness=0.5));
      connect(H2_aq.outlet, gasVolatility.inlet) annotation (Line(
          points={{44,-66},{44,-36},{40,-36}},
          color={158,66,200},
          thickness=0.5));
      connect(gasVolatility.outlet, H2.inlet) annotation (Line(
          points={{60,-36},{66,-36},{66,42},{60,42}},
          color={200,66,175},
          thickness=0.5));
      connect(gases.solution, gasVolatility.solution) annotation (Line(points={{39.6,22.36},{39.6,18},{60.2,18},{60.2,-40.2}}, color={127,127,0}));
      connect(solution1.solution, gasSolubility.solution) annotation (Line(points={{38.8,-95.34},{38.8,-100},{28,-100},{28,-96},{-2,-96},{-2,-102},{-42,-102},{
              -42,-26},{-12,-26},{-12,-16.2},{-15.8,-16.2}}, color={127,127,0}));
      annotation (
      experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
    end H2O_ElectrochemicalCell;

    package ClimateChange
      class References "References"
        extends Modelica.Icons.References;

        annotation (DocumentationClass=true,Documentation(info="<html>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"0\">
<tr>
<td>[Nelabhotla2019]</td>
<td>Anirudh Bhanu Teja Nelabhotla, Rune Bakke, Carlos Dinamarca,
        \"Performance Analysis of Biocathode in Bioelectrochemical CO2 Reduction\"
         Catalysts, 9, 683, 2019,
        <a href=\"https://doi.org/10.3390/catal9080683\">doi:10.3390/catal9080683</a>.
</td>
</tr>
</table>
</html>"));
      end References;

      model AcetoclasticMethanogenesis
        extends Modelica.Icons.Example;
        Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
        Chemical.Boundaries.Substance CH3COOH(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.AceticAcid_aqueous(),
          mass_start=0.001) "Acetic acid" annotation (Placement(transformation(extent={{-72,30},{-52,50}})));
        Chemical.Boundaries.Substance CH4(
          use_mass_start=false,
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{26,48},{46,68}})));
        Chemical.Boundaries.Substance CO2(
          use_mass_start=false,
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          mass_start=0.001) annotation (Placement(transformation(extent={{56,10},{76,30}})));
        Chemical.Boundaries.Substance Water(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{-14,-66},{6,-46}})));
        Chemical.Processes.Reaction reaction(
          productsSubstanceData={Chemical.SubstancesOld.Methan_aqueous(),Chemical.SubstancesOld.CarbonDioxide_aqueous()},
          nS=1,
          nP=2) "Acetoclastic (heterotrophic) methanogenesis" annotation (Placement(transformation(extent={{-18,32},{2,52}})));
        inner DropOfCommons dropOfCommons(L=10) annotation (Placement(transformation(extent={{-40,-26},{-20,-6}})));
      equation
        connect(CH3COOH.solution, solution.solution) annotation (Line(points={{
                -68,30},{-68,-88},{60,-88},{60,-98}}, color={127,127,0}));
        connect(Water.solution, solution.solution) annotation (Line(points={{-10,
                -66},{-10,-88},{60,-88},{60,-98}}, color={127,127,0}));
        connect(CO2.solution, solution.solution)
          annotation (Line(points={{60,10},{60,-98}}, color={127,127,0}));
        connect(CH4.solution, solution.solution) annotation (Line(points={{30,48},
                {32,48},{32,-88},{60,-88},{60,-98}}, color={127,127,0}));
        connect(CH3COOH.outlet, reaction.substrates[1]) annotation (Line(
            points={{-52,40},{-36,40},{-36,42},{-18,42}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[1], CH4.inlet) annotation (Line(
            points={{2,41.75},{16,41.75},{16,58},{26,58}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[2], CO2.inlet)
          annotation (Line(
            points={{2,42.25},{2,38},{16,38},{16,20},{56,20}},
            color={158,66,200},
            thickness=0.5));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=10, __Dymola_Algorithm="Dassl"),
          __Dymola_experimentSetupOutput);
      end AcetoclasticMethanogenesis;

      model HydrogenotrophicMethanogenesis
        extends Modelica.Icons.Example;
        Solution solution annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
        Chemical.Boundaries.Substance CH4(
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{32,72},{52,92}})));
        Chemical.Boundaries.Substance CO2(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
          mass_start=0.001) annotation (Placement(transformation(extent={{-68,40},{-48,60}})));
        Chemical.Boundaries.Substance H2O(
          use_mass_start=false,
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          mass_start=1) annotation (Placement(transformation(extent={{36,36},{56,56}})));
        Chemical.Processes.Reaction reaction(
          s={4,1},
          p={1,2},
          productsSubstanceData={Chemical.SubstancesOld.Methan_aqueous(),Chemical.SubstancesOld.Water_liquid()},
          nS=2,
          nP=2) "Hydrogenotrophic (autotrophic) methanogenesis" annotation (Placement(transformation(extent={{-18,50},{2,70}})));
        Chemical.Boundaries.Substance H2(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.Hydrogen_aqueous(),
          mass_start=0.001) annotation (Placement(transformation(extent={{-74,66},{-54,86}})));
        inner DropOfCommons dropOfCommons(L=1e12) annotation (Placement(transformation(extent={{-64,-56},{-44,-36}})));
      equation
        connect(H2O.solution, solution.solution) annotation (Line(points={{40,36},
                {40,22},{60,22},{60,-98}}, color={127,127,0}));
        connect(CO2.solution, solution.solution) annotation (Line(points={{-64,40},
                {-64,22},{60,22},{60,-98}}, color={127,127,0}));
        connect(CH4.solution, solution.solution) annotation (Line(points={{36,72},{36,22},{60,22},{60,-98}},
                                           color={127,127,0}));
        connect(H2.solution, solution.solution) annotation (Line(points={{-70,66},
                {-70,22},{60,22},{60,-98}}, color={127,127,0}));
        connect(H2.outlet, reaction.substrates[1])
          annotation (Line(
            points={{-54,76},{-38,76},{-38,60},{-18,60},{-18,59.75}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.outlet, reaction.substrates[2])
          annotation (Line(
            points={{-48,50},{-34,50},{-34,56},{-18,56},{-18,60.25}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[1], CH4.inlet) annotation (Line(
            points={{2,59.75},{2,60},{32,60},{32,82}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[2], H2O.inlet) annotation (Line(
            points={{2,60.25},{2,58},{36,58},{36,46}},
            color={158,66,200},
            thickness=0.5));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StopTime=2352,
            Tolerance=1e-08,
            __Dymola_Algorithm="Dassl"),
          __Dymola_experimentSetupOutput);
      end HydrogenotrophicMethanogenesis;

      model MethanElectrosynthesis
        "Direct electron transfer (electrosynthesis reaction-bioelectrochemical methane)"
       extends Modelica.Icons.Example;

        Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
        Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

        Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-96},{38,6}})));

        Chemical.Boundaries.ExternalIdealGas CH4(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Methan_gas(),
          PartialPressure=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={34,44})));
        Chemical.Boundaries.Substance H(
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-16,-28},{4,-8}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-20,62},{0,82}})));
        Chemical.Processes.Reaction electrodeReaction(
          s={1,8,8},
          p={1,2},
          nS=3,
          nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={52,6})));
        Chemical.Processes.Reaction electrodeReaction1(
          s={2,8,8},
          p={4},
          nS=3,
          nP=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-40,0})));

        Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{46,52},{66,72}})));
        Chemical.Boundaries.Substance H2O(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{-6,-80},{14,-60}})));
        Chemical.Boundaries.ExternalIdealGas O2_(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
        Chemical.Boundaries.ExternalIdealGas CO2(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{20,36},{0,56}})));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
          annotation (Placement(transformation(extent={{-80,78},{-60,98}})));
        Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
          annotation (Placement(transformation(extent={{12,78},{32,98}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
          annotation (Placement(transformation(extent={{76,78},{96,98}})));
        Chemical.Topology.SplitterT1 splitterT1
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={4,14})));
        Chemical.Topology.JunctionT1 junctionT1
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-10,-46})));
        Solution gases(BasePressure=100000) annotation (Placement(transformation(extent={{-36,22},{50,64}})));
      equation
        connect(H.solution, solution1.solution) annotation (Line(points={{-12,-28},{-12,-32},{24.4,-32},{24.4,-94.98}},
                                           color={127,127,0}));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-74,32},{-74,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
          color={127,127,0}));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{84,-26},{84,-49},{89.2,-49}},
          color={127,127,0}));
        connect(voltageSensor.p, electrone.pin) annotation (Line(
            points={{-20,72},{-88,72},{-88,51.8},{-68,51.8}},
            color={0,0,255}));
        connect(voltageSensor.n, electrone1.pin) annotation (Line(
            points={{0,72},{98,72},{98,-6.2},{78,-6.2}},
            color={0,0,255}));
        connect(electrone1.pin, ground.p) annotation (Line(
            points={{78,-6.2},{98,-6.2},{98,72},{56,72}},
            color={0,0,255}));
        connect(H2O.solution, solution1.solution) annotation (Line(points={{-2,-80},
                {-2,-94.98},{24.4,-94.98}}, color={127,127,0}));
        connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-68,51.8},{-88,51.8},{-88,88},{-80,88}},
                                                   color={0,0,255}));
        connect(constantVoltage.n, resistor.p)
          annotation (Line(points={{-60,88},{12,88}},   color={0,0,255}));
        connect(powerSensor.nv, resistor.n) annotation (Line(points={{86,78},{74,
                78},{74,88},{32,88}},   color={0,0,255}));
        connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{86,98},{-80,98},{-80,88}},
                                                     color={0,0,255}));
        connect(resistor.n, powerSensor.pc) annotation (Line(points={{32,88},{76,88}},
                                          color={0,0,255}));
        connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{96,88},{98,88},{98,-6.2},{78,-6.2}},
                                               color={0,0,255}));
        connect(O2_.outlet, electrodeReaction1.substrates[1])
          annotation (Line(
            points={{-26,42},{-39.6667,42},{-39.6667,10}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletA, electrodeReaction1.substrates[2])
          annotation (Line(
            points={{-6,14},{-40,14},{-40,10}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone.outlet, electrodeReaction1.substrates[3])
          annotation (Line(
            points={{-58,42},{-40.3333,42},{-40.3333,10}},
            color={158,66,200},
            thickness=0.5));
        connect(H.outlet, splitterT1.inlet) annotation (Line(
            points={{4,-18},{4,4}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.outlet, electrodeReaction.substrates[1]) annotation (Line(
            points={{0,46},{0,18},{32,18},{32,-12},{52,-12},{52,-18},{51.6667,-18},{51.6667,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletB, electrodeReaction.substrates[2])
          annotation (Line(
            points={{14,14},{28,14},{28,-20},{52,-20},{52,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone1.outlet, electrodeReaction.substrates[3])
          annotation (Line(
            points={{68,-16},{60,-16},{60,-18},{52.3333,-18},{52.3333,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction.products[1], CH4.inlet)
          annotation (Line(
            points={{51.75,16},{51.75,44},{44,44},{44,44}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1.outlet, H2O.inlet) annotation (Line(
            points={{-10,-56},{-10,-70},{-6,-70}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[1], junctionT1.inletA)
          annotation (Line(
            points={{-40,-10},{-40,-46},{-20,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1.inletB, electrodeReaction.products[2])
          annotation (Line(
            points={{0,-46},{36,-46},{36,24},{52.25,24},{52.25,16}},
            color={158,66,200},
            thickness=0.5));
        connect(gases.solution, O2_.solution) annotation (Line(points={{32.8,22.42},{11.4,22.42},{11.4,32},{-10,32}}, color={127,127,0}));
        connect(gases.solution, CH4.solution) annotation (Line(points={{32.8,22.42},{32.8,28.21},{40,28.21},{40,34}}, color={127,127,0}));
        connect(CO2.solution, gases.solution) annotation (Line(points={{16,36},{26,36},{26,22.42},{32.8,22.42}}, color={127,127,0}));
        annotation (
        experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end MethanElectrosynthesis;

      model ElectrolysisMethanation
       extends Modelica.Icons.Example;

        Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
        Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,26},{94,72}})));

        Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

        Chemical.Boundaries.Substance H(
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Processes.Reaction electrodeReaction(
          s={2,2},
          nS=2,
          nP=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={42,34})));
        Chemical.Processes.Reaction electrodeReaction1(
          s={2,8,8},
          p={4},
          nS=3,
          nP=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-36,28})));

        Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{22,54},{42,74}})));
        Chemical.Boundaries.Substance H2O(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{68,-28},{48,-8}})));
        Chemical.Boundaries.ExternalIdealGas O2_(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{14,40},{-6,60}})));
        Chemical.Boundaries.Substance CH4(
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.Methan_aqueous(),
          mass_start=0.001) "Methan" annotation (Placement(transformation(extent={{54,-64},{74,-44}})));
        Chemical.Boundaries.Substance CO2(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
          mass_start=0.1) annotation (Placement(transformation(extent={{-72,-76},{-52,-56}})));
        Chemical.Processes.Reaction reaction(
          s={4,1},
          p={1,2},
          nP=2,
          nS=2) "Hydrogenotrophic (autotrophic) methanogenesis" annotation (Placement(transformation(extent={{-2,-54},{18,-34}})));
        Chemical.Boundaries.Substance H2(
          useInlet=true,
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.Hydrogen_aqueous(),
          mass_start=5e-11) annotation (Placement(transformation(extent={{-76,-44},{-56,-24}})));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
          annotation (Placement(transformation(extent={{-72,84},{-52,104}})));
        Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
          annotation (Placement(transformation(extent={{20,84},{40,104}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
          annotation (Placement(transformation(extent={{112,80},{132,100}})));
        Solution gas(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas                                          "Ideal Gas", BasePressure=100000)
          annotation (Placement(transformation(extent={{-24,26},{22,68}})));
        Chemical.Topology.SplitterT2 splitterT2 annotation (Placement(transformation(extent={{-52,-16},{-32,4}})));
        Chemical.Topology.JunctionT1 junctionT1 annotation (Placement(transformation(extent={{72,-28},{92,-8}})));
      equation
        connect(H.solution, solution1.solution) annotation (Line(points={{-80,-14},{-80,
                -88},{55.6,-88},{55.6,-88.96}},
                                           color={127,127,0}));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-74,44},{-54,44},{-54,26.46}},
          color={127,127,0}));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{84,38},{84,26.46},{86.8,26.46}},
          color={127,127,0}));
        connect(voltageSensor.p, electrone.pin) annotation (Line(
            points={{-6,74},{-96,74},{-96,63.8},{-68,63.8}},
            color={0,0,255}));
        connect(voltageSensor.n, electrone1.pin) annotation (Line(
            points={{14,74},{92,74},{92,57.8},{78,57.8}},
            color={0,0,255}));
        connect(electrone1.pin, ground.p) annotation (Line(
            points={{78,57.8},{92,57.8},{92,74},{32,74}},
            color={0,0,255}));
        connect(H2O.solution, solution1.solution) annotation (Line(points={{64,-28},{64,-88.96},{55.6,-88.96}},
                                            color={127,127,0}));
        connect(CO2.solution, solution1.solution) annotation (Line(points={{-68,-76},{
                -80,-76},{-80,-88.96},{55.6,-88.96}}, color={127,127,0}));
        connect(H2.solution, solution1.solution) annotation (Line(points={{-72,-44},{-80,
                -44},{-80,-88},{56,-88},{56,-88.96},{55.6,-88.96}}, color={127,127,0}));
        connect(CH4.solution, solution1.solution) annotation (Line(points={{58,-64},{74,-64},{74,-88.96},{55.6,-88.96}},
                                                 color={127,127,0}));
        connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-68,63.8},{-96,63.8},{-96,94},{-72,94}},
                                                 color={0,0,255}));
        connect(constantVoltage.n, resistor.p)
          annotation (Line(points={{-52,94},{20,94}}, color={0,0,255}));
        connect(powerSensor.nv, resistor.n) annotation (Line(points={{122,80},{82,
                80},{82,94},{40,94}}, color={0,0,255}));
        connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{122,
                100},{-74,100},{-74,94},{-72,94}}, color={0,0,255}));
        connect(resistor.n, powerSensor.pc) annotation (Line(points={{40,94},{76,
                94},{76,90},{112,90}}, color={0,0,255}));
        connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{132,90},{142,90},{142,57.8},{78,57.8}},
                                            color={0,0,255}));
        connect(gas.solution, O2_.solution) annotation (Line(points={{12.8,26.42},{10,26.42},{10,40}},        color={127,127,0}));
        connect(O2_.outlet, electrodeReaction1.substrates[1])
          annotation (Line(
            points={{-6,50},{-35.6667,50},{-35.6667,38}},
            color={158,66,200},
            thickness=0.5));
        connect(H.outlet, splitterT2.inletProcess) annotation (Line(
            points={{-64,-4},{-64,-6},{-52,-6}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT2.outletSubstance, electrodeReaction1.substrates[2])
          annotation (Line(
            points={{-42,4},{-42,38},{-36,38}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone.outlet, electrodeReaction1.substrates[3])
          annotation (Line(
            points={{-58,54},{-48,54},{-48,58},{-36.3333,58},{-36.3333,38}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT2.outletSubstance1, electrodeReaction.substrates[1])
          annotation (Line(
            points={{-32,-6},{41.75,-6},{41.75,24}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone1.outlet, electrodeReaction.substrates[2])
          annotation (Line(
            points={{68,48},{58,48},{58,2},{42.25,2},{42.25,24}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction.products[1], H2.inlet)
          annotation (Line(
            points={{42,44},{40,44},{40,56},{26,56},{26,-20},{-108,-20},{-108,-34},{-76,-34}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[1], CH4.inlet)
          annotation (Line(
            points={{18,-44.25},{38,-44.25},{38,-54},{54,-54}},
            color={158,66,200},
            thickness=0.5));
        connect(H2O.inlet, junctionT1.outlet) annotation (Line(
            points={{68,-18},{72,-18}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[1], junctionT1.inletA)
          annotation (Line(
            points={{-36,18},{-36,6},{82,6},{82,-8}},
            color={158,66,200},
            thickness=0.5));
        connect(reaction.products[2], junctionT1.inletB)
          annotation (Line(
            points={{18,-43.75},{18,-38},{82,-38},{82,-28}},
            color={158,66,200},
            thickness=0.5));
        connect(H2.outlet, reaction.substrates[1])
          annotation (Line(
            points={{-56,-34},{-44,-34},{-44,-40},{-2,-40},{-2,-44.25}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.outlet, reaction.substrates[2])
          annotation (Line(
            points={{-52,-66},{-28,-66},{-28,-54},{-2,-54},{-2,-43.75}},
            color={158,66,200},
            thickness=0.5));
        annotation (
        experiment(StopTime=0.0001, __Dymola_Algorithm="Dassl"),
                                     Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end ElectrolysisMethanation;

      model ElectrochemicalAcetateProduction
        "Direct electron transfer (electrochemical acetate production)"
       extends Modelica.Icons.Example;

        Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-86,26},{-46,72}})));
        Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{58,26},{94,72}})));

        Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-90,-90},{92,14}})));

        Chemical.Boundaries.Substance H(
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-4) annotation (Placement(transformation(extent={{-84,-14},{-64,6}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-6,64},{14,84}})));
        Chemical.Processes.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-36,28})));

        Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,44},{-58,64}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,38},{68,58}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{22,54},{42,74}})));
        Chemical.Boundaries.Substance H2O(substanceData=Chemical.SubstancesOld.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{78,-28},{58,-8}})));
        Chemical.Boundaries.ExternalIdealGas O2_(
          substanceData=Chemical.SubstancesOld.Oxygen_gas(),
          PartialPressure=100000,
          TotalPressure=100000) annotation (Placement(transformation(extent={{0,36},{-20,56}})));
        Chemical.Boundaries.Substance AcAc(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.AceticAcid_aqueous(),
          mass_start=0.001) "Acetic Acid" annotation (Placement(transformation(extent={{60,-60},{40,-40}})));
        Chemical.Boundaries.Substance CO2(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.SubstancesOld.CarbonDioxide_aqueous(),
          mass_start=1) annotation (Placement(transformation(extent={{-10,-34},{10,-14}})));
        Chemical.Processes.Reaction reaction(
          s={2,8,8},
          p={1,2},
          nS=3,
          nP=2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={50,30})));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
          annotation (Placement(transformation(extent={{-54,76},{-34,96}})));
        Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
          annotation (Placement(transformation(extent={{38,76},{58,96}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
          annotation (Placement(transformation(extent={{130,72},{150,92}})));
      equation
        connect(H.solution, solution1.solution) annotation (Line(points={{-80,-14},{-80,
                -88},{55.6,-88},{55.6,-88.96}},
                                           color={127,127,0}));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-74,44},{-54,44},{-54,26.46}},
          color={127,127,0}));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{84,38},{84,26.46},{86.8,26.46}},
          color={127,127,0}));
        connect(voltageSensor.p, electrone.pin) annotation (Line(
            points={{-6,74},{-96,74},{-96,63.8},{-68,63.8}},
            color={0,0,255}));
        connect(voltageSensor.n, electrone1.pin) annotation (Line(
            points={{14,74},{92,74},{92,57.8},{78,57.8}},
            color={0,0,255}));
        connect(electrone1.pin, ground.p) annotation (Line(
            points={{78,57.8},{92,57.8},{92,74},{32,74}},
            color={0,0,255}));
        connect(H2O.solution, solution1.solution) annotation (Line(points={{74,-28},{74,
                -88.96},{55.6,-88.96}},     color={127,127,0}));
        connect(electrodeReaction1.substrates[1], H2O.port_a) annotation (Line(
              points={{-36,18},{-36,2},{38,2},{38,-18},{58,-18}},
                                                     color={158,66,200}));
        connect(O2_.port_a, electrodeReaction1.products[1]) annotation (Line(points={{-20,46},{-38,46},{-38,38},{-37.3333,38}},
                                                          color={158,66,200}));
        connect(H.port_a, electrodeReaction1.products[2]) annotation (Line(points={{-64,-4},
                {-28,-4},{-28,58},{-36,58},{-36,38}},   color={158,66,200}));
        connect(electrone.port_a, electrodeReaction1.products[3]) annotation (Line(
              points={{-58,54},{-42,54},{-42,38},{-34.6667,38}}, color={158,66,200}));
        connect(CO2.solution, solution1.solution) annotation (Line(points={{-6,-34},
                {-80,-34},{-80,-88.96},{55.6,-88.96}},color={127,127,0}));
        connect(AcAc.solution, solution1.solution) annotation (Line(points={{56,
                -60},{58,-60},{58,-88.96},{55.6,-88.96}}, color={127,127,0}));
        connect(CO2.port_a, reaction.substrates[1]) annotation (Line(points={{10,-24},{20,-24},{20,52},{48.6667,52},{48.6667,40}},
                                                                  color={158,66,
                200}));
        connect(H.port_a, reaction.substrates[2]) annotation (Line(points={{-64,
                -4},{14,-4},{14,56},{50,56},{50,40}}, color={158,66,200}));
        connect(electrone1.port_a, reaction.substrates[3]) annotation (Line(
              points={{68,48},{46,48},{46,40},{51.3333,40}}, color={158,66,200}));
        connect(AcAc.port_a, reaction.products[1]) annotation (Line(points={{40,-50},{24,-50},{24,8},{49,8},{49,20}},
                                                      color={158,66,200}));
        connect(H2O.port_a, reaction.products[2]) annotation (Line(points={{58,-18},{58,-19},{51,-19},{51,20}},
                                                 color={158,66,200}));
        connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-68,63.8},{-96,63.8},{-96,86},{-54,86}},
                                                 color={0,0,255}));
        connect(constantVoltage.n, resistor.p)
          annotation (Line(points={{-34,86},{38,86}}, color={0,0,255}));
        connect(powerSensor.nv, resistor.n) annotation (Line(points={{140,72},{
                100,72},{100,86},{58,86}}, color={0,0,255}));
        connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{140,
                92},{-56,92},{-56,86},{-54,86}}, color={0,0,255}));
        connect(resistor.n, powerSensor.pc) annotation (Line(points={{58,86},{94,
                86},{94,82},{130,82}}, color={0,0,255}));
        connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{150,82},{162,82},{162,57.8},{78,57.8}},
                                            color={0,0,255}));
        annotation (
        experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end ElectrochemicalAcetateProduction;

      model MethanElectrosynthesis_ "Direct electron transfer (electrosynthesis reaction-bioelectrochemical methane)"
       extends Modelica.Icons.Example;

        Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
        Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

        Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-96},{38,6}})));

        Chemical.Boundaries.ExternalIdealGas CH4(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Methan_gas(),
          PartialPressure=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={34,44})));
        Chemical.Boundaries.Substance H(
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-16,-28},{4,-8}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-20,62},{0,82}})));
        Chemical.Processes.Reaction electrodeReaction(
          s={1,8,8},
          p={1,2},
          nS=3,
          nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={52,6})));
        Chemical.Processes.Reaction electrodeReaction1(
          s={2,8,8},
          p={4},
          nS=3,
          nP=1) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-40,0})));

        Chemical.Boundaries.ElectronSource electrone annotation (Placement(transformation(extent={{-78,32},{-58,52}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{46,52},{66,72}})));
        Chemical.Boundaries.Substance H2O(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{-6,-80},{14,-60}})));
        Chemical.Boundaries.ExternalIdealGas O2_(substanceData=Chemical.SubstancesOld.Oxygen_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{-6,32},{-26,52}})));
        Chemical.Boundaries.ExternalIdealGas CO2(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{20,36},{0,56}})));
        Chemical.Topology.SplitterT1 splitterT1
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={4,14})));
        Chemical.Topology.JunctionT1 junctionT1
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-10,-46})));
        Solution gases(BasePressure=100000) annotation (Placement(transformation(extent={{-36,22},{50,64}})));
      equation
        connect(H.solution, solution1.solution) annotation (Line(points={{-12,-28},{-12,-32},{24.4,-32},{24.4,-94.98}},
                                           color={127,127,0}));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-74,32},{-74,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
          color={127,127,0}));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{84,-26},{84,-49},{89.2,-49}},
          color={127,127,0}));
        connect(voltageSensor.p, electrone.pin) annotation (Line(
            points={{-20,72},{-88,72},{-88,51.8},{-68,51.8}},
            color={0,0,255}));
        connect(voltageSensor.n, electrone1.pin) annotation (Line(
            points={{0,72},{98,72},{98,-6.2},{78,-6.2}},
            color={0,0,255}));
        connect(electrone1.pin, ground.p) annotation (Line(
            points={{78,-6.2},{98,-6.2},{98,72},{56,72}},
            color={0,0,255}));
        connect(H2O.solution, solution1.solution) annotation (Line(points={{-2,-80},
                {-2,-94.98},{24.4,-94.98}}, color={127,127,0}));
        connect(O2_.outlet, electrodeReaction1.substrates[1])
          annotation (Line(
            points={{-26,42},{-39.6667,42},{-39.6667,10}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletA, electrodeReaction1.substrates[2])
          annotation (Line(
            points={{-6,14},{-40,14},{-40,10}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone.outlet, electrodeReaction1.substrates[3])
          annotation (Line(
            points={{-58,42},{-40.3333,42},{-40.3333,10}},
            color={158,66,200},
            thickness=0.5));
        connect(H.outlet, splitterT1.inlet) annotation (Line(
            points={{4,-18},{4,4}},
            color={158,66,200},
            thickness=0.5));
        connect(CO2.outlet, electrodeReaction.substrates[1]) annotation (Line(
            points={{0,46},{0,18},{32,18},{32,-12},{52,-12},{52,-18},{51.6667,-18},{51.6667,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(splitterT1.outletB, electrodeReaction.substrates[2])
          annotation (Line(
            points={{14,14},{28,14},{28,-20},{52,-20},{52,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone1.outlet, electrodeReaction.substrates[3])
          annotation (Line(
            points={{68,-16},{60,-16},{60,-18},{52.3333,-18},{52.3333,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction.products[1], CH4.inlet)
          annotation (Line(
            points={{51.75,16},{51.75,44},{44,44},{44,44}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1.outlet, H2O.inlet) annotation (Line(
            points={{-10,-56},{-10,-70},{-6,-70}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[1], junctionT1.inletA)
          annotation (Line(
            points={{-40,-10},{-40,-46},{-20,-46}},
            color={158,66,200},
            thickness=0.5));
        connect(junctionT1.inletB, electrodeReaction.products[2])
          annotation (Line(
            points={{0,-46},{36,-46},{36,24},{52.25,24},{52.25,16}},
            color={158,66,200},
            thickness=0.5));
        connect(gases.solution, O2_.solution) annotation (Line(points={{32.8,22.42},{11.4,22.42},{11.4,32},{-10,32}}, color={127,127,0}));
        connect(gases.solution, CH4.solution) annotation (Line(points={{32.8,22.42},{32.8,28.21},{40,28.21},{40,34}}, color={127,127,0}));
        connect(CO2.solution, gases.solution) annotation (Line(points={{16,36},{26,36},{26,22.42},{32.8,22.42}}, color={127,127,0}));
        annotation (
        experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end MethanElectrosynthesis_;

      model MethanElectrosynthesis2 "Direct electron transfer (electrosynthesis reaction-bioelectrochemical methane)"
       extends Modelica.Icons.Example;

        Chemical.Solution cathode(ElectricGround=false) annotation (Placement(transformation(extent={{-88,-44},{-46,72}})));
        Chemical.Solution anode(ElectricGround=false) annotation (Placement(transformation(extent={{62,-50},{96,50}})));

        Chemical.Solution solution1(ElectricGround=false) annotation (Placement(transformation(extent={{-30,-96},{38,6}})));

        Chemical.Boundaries.ExternalIdealGas CH4(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Methan_gas(),
          PartialPressure=100000) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={34,44})));
        Chemical.Boundaries.Substance H(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Proton_aqueous(),
          use_mass_start=false,
          amountOfSubstance_start=1e-7) annotation (Placement(transformation(extent={{-16,-28},{4,-8}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(extent={{-20,62},{0,82}})));
        Chemical.Processes.Reaction electrodeReaction(
          s={1,8,8},
          p={1,2},
          nS=3,
          nP=2) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={52,6})));
        Chemical.Processes.Reaction electrodeReaction1(
          s={4},
          p={2,8,8},
          nS=1,
          nP=3) annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=270,
              origin={-40,0})));

        Chemical.Boundaries.ElectronSink electrone annotation (Placement(transformation(extent={{-58,32},{-78,52}})));
                                   //(substanceData=Chemical.Examples.Substances.Electrone_solid())
        Chemical.Boundaries.ElectronSource electrone1 annotation (Placement(transformation(extent={{88,-26},{68,-6}})));
                                    //(substanceData=Chemical.Examples.Substances.Electrone_solid())
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{-112,62},{-92,82}})));
        Chemical.Boundaries.Substance H2O(
          useInlet=true,
          substanceData=Chemical.SubstancesOld.Water_liquid(),
          mass_start=1) annotation (Placement(transformation(extent={{14,-80},{-6,-60}})));
        Chemical.Boundaries.Substance O2_(useInlet=true) annotation (Placement(transformation(extent={{-26,32},{-6,52}})));
        Chemical.Boundaries.ExternalIdealGas CO2(substanceData=Chemical.SubstancesOld.CarbonDioxide_gas(), PartialPressure=100000)
          annotation (Placement(transformation(extent={{20,36},{0,56}})));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=12)
          annotation (Placement(transformation(extent={{-80,78},{-60,98}})));
        Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.01)
          annotation (Placement(transformation(extent={{12,78},{32,98}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
          annotation (Placement(transformation(extent={{76,78},{96,98}})));
        Solution gases(redeclare package stateOfMatter =
              Chemical.Interfaces.IdealGas                                            "Ideal Gas", BasePressure=100000)
          annotation (Placement(transformation(extent={{-36,22},{50,64}})));
      equation
        connect(H.solution, solution1.solution) annotation (Line(points={{-12,-28},{-12,-32},{24.4,-32},{24.4,-94.98}},
                                           color={127,127,0}));
      connect(electrone.solution, cathode.solution) annotation (Line(
          points={{-62,32},{-62,-18},{-64,-18},{-64,-42.84},{-54.4,-42.84}},
          color={127,127,0}));
      connect(electrone1.solution, anode.solution) annotation (Line(
          points={{84,-26},{84,-49},{89.2,-49}},
          color={127,127,0}));
        connect(voltageSensor.p, electrone.pin) annotation (Line(
            points={{-20,72},{-88,72},{-88,51.8},{-68,51.8}},
            color={0,0,255}));
        connect(voltageSensor.n, electrone1.pin) annotation (Line(
            points={{0,72},{98,72},{98,-6.2},{78,-6.2}},
            color={0,0,255}));
        connect(H2O.solution, solution1.solution) annotation (Line(points={{10,-80},{10,-94.98},{24.4,-94.98}},
                                            color={127,127,0}));
        connect(electrone.pin, constantVoltage.p) annotation (Line(points={{-68,51.8},{-88,51.8},{-88,88},{-80,88}},
                                                   color={0,0,255}));
        connect(constantVoltage.n, resistor.p)
          annotation (Line(points={{-60,88},{12,88}},   color={0,0,255}));
        connect(powerSensor.nv, resistor.n) annotation (Line(points={{86,78},{74,
                78},{74,88},{32,88}},   color={0,0,255}));
        connect(powerSensor.pv, constantVoltage.p) annotation (Line(points={{86,98},{-80,98},{-80,88}},
                                                     color={0,0,255}));
        connect(resistor.n, powerSensor.pc) annotation (Line(points={{32,88},{76,88}},
                                          color={0,0,255}));
        connect(powerSensor.nc, electrone1.pin) annotation (Line(points={{96,88},{98,88},{98,-6.2},{78,-6.2}},
                                               color={0,0,255}));
        connect(CO2.outlet, electrodeReaction.substrates[1]) annotation (Line(
            points={{0,46},{0,18},{32,18},{32,-12},{52,-12},{52,-18},{51.6667,-18},{51.6667,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrone1.outlet, electrodeReaction.substrates[2])
          annotation (Line(
            points={{68,-16},{60,-16},{60,-18},{52,-18},{52,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction.products[1], CH4.inlet)
          annotation (Line(
            points={{51.75,16},{51.75,44},{44,44},{44,44}},
            color={158,66,200},
            thickness=0.5));
        connect(gases.solution, O2_.solution) annotation (Line(points={{32.8,22.42},{11.4,22.42},{11.4,32},{-22,32}}, color={127,127,0}));
        connect(gases.solution, CH4.solution) annotation (Line(points={{32.8,22.42},{32.8,28.21},{40,28.21},{40,34}}, color={127,127,0}));
        connect(CO2.solution, gases.solution) annotation (Line(points={{16,36},{26,36},{26,22.42},{32.8,22.42}}, color={127,127,0}));
        connect(H2O.inlet, electrodeReaction.products[2])
          annotation (Line(
            points={{14,-70},{36,-70},{36,24},{52.25,24},{52.25,16}},
            color={158,66,200},
            thickness=0.5));
        connect(H2O.outlet, electrodeReaction1.substrates[1])
          annotation (Line(
            points={{-6,-70},{-38,-70},{-38,-64},{-40,-64},{-40,-10}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[1], O2_.inlet)
          annotation (Line(
            points={{-39.6667,10},{-39.6667,46},{-26,46},{-26,42}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[2], electrone.inlet)
          annotation (Line(
            points={{-40,10},{-40,42},{-78,42}},
            color={158,66,200},
            thickness=0.5));
        connect(electrodeReaction1.products[3], H.inlet)
          annotation (Line(
            points={{-40.3333,10},{-28,10},{-28,-18},{-16,-18}},
            color={158,66,200},
            thickness=0.5));
        connect(H.outlet, electrodeReaction.substrates[3])
          annotation (Line(
            points={{4,-18},{52.3333,-18},{52.3333,-4}},
            color={158,66,200},
            thickness=0.5));
        connect(ground.p, constantVoltage.p) annotation (Line(points={{-102,82},{-102,88},{-80,88}}, color={0,0,255}));
        annotation (
        experiment(StopTime=1),      Documentation(revisions="<html>
<p><i>2015-2018</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"));
      end MethanElectrosynthesis2;
    end ClimateChange;

    model PKPD
      Chemical.Boundaries.TerminalInflow substanceInflow(SubstanceFlow=2) annotation (Placement(transformation(extent={{-82,22},{-62,42}})));
      Chemical.Boundaries.Substance pumped(useInlet=true, useSolution=false) annotation (Placement(transformation(extent={{-28,22},{-8,42}})));
      Chemical.Boundaries.TerminalOutflow outflow(SubstanceFlow=1) annotation (Placement(transformation(extent={{54,22},{74,42}})));
      Chemical.Boundaries.TerminalInflow substanceInflow2(SubstanceFlow=2) annotation (Placement(transformation(extent={{-82,-10},{-62,10}})));
      Chemical.Boundaries.Substance clearanced(useInlet=true, useSolution=false) annotation (Placement(transformation(extent={{-28,-10},{-8,10}})));
      Chemical.Boundaries.Clearance clearance(Clearance(displayUnit="l/s") = 0.002) annotation (Placement(transformation(extent={{54,-10},{74,10}})));
      Chemical.Boundaries.Degradation degradation(HalfTime(displayUnit="min") = 60) annotation (Placement(transformation(extent={{54,-48},{74,-28}})));
      Chemical.Boundaries.TerminalInflow substanceInflow3(SubstanceFlow=2) annotation (Placement(transformation(extent={{-82,-48},{-62,-28}})));
      Chemical.Boundaries.Substance degraded(useInlet=true, useSolution=false) annotation (Placement(transformation(extent={{-28,-48},{-8,-28}})));
    equation
      connect(substanceInflow.outlet, pumped.inlet) annotation (Line(
          points={{-62,32},{-28,32}},
          color={158,66,200},
          thickness=0.5));
      connect(pumped.outlet, outflow.inlet) annotation (Line(
          points={{-8,32},{54,32}},
          color={158,66,200},
          thickness=0.5));
      connect(substanceInflow2.outlet, clearanced.inlet) annotation (Line(
          points={{-62,0},{-28,0}},
          color={158,66,200},
          thickness=0.5));
      connect(substanceInflow3.outlet, degraded.inlet) annotation (Line(
          points={{-62,-38},{-28,-38}},
          color={158,66,200},
          thickness=0.5));
      connect(clearanced.outlet, clearance.inlet) annotation (Line(
          points={{-8,0},{54,0}},
          color={158,66,200},
          thickness=0.5));
      connect(degraded.outlet, degradation.inlet) annotation (Line(
          points={{-8,-38},{54,-38}},
          color={158,66,200},
          thickness=0.5));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Dassl"));
    end PKPD;
  end Examples;
end Onedirectional;
