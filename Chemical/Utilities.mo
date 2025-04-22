within Chemical;
package Utilities
  package Units
    type Inertance =  Real (quantity = "Inertance",unit = "s2/mol", start=1e-5, nominal=1e-5, min=0) "Inertance of electro-chemical process";
    type URT =  Real (quantity = "URT",unit = "1", nominal=1)  "Electro-chemical potential divided by gas constant and temperature";
    type MolarFlowAcceleration = Real(quantity="MolarFlowAcceleration", unit="mol/s2") "Acceleration Unit for a MolarFlow"
      annotation (Documentation(info="<html>
    <p>Acceleration Unit for a MolarFlow</p>
    </html>"));
  end Units;

  package Types
    type InitializationMethods = enumeration(
      none
        "No initialization",
      steadyState
        "Steady state initialization (derivatives of states are zero)",
      state
        "Initialization with initial states",
      derivative
        "Initialization with initial derivatives of states") "Choices for initialization of a state."
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(
          coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>
Choices for initialization of a&nbsp;state.
</p>
</html>"));
    type InitializationSubstance
                               = enumeration(
      none
        "No initialization",
      steadyStateForwards
        "Steady state initialization in forward direction (chemical potential from rear connector)",
      steadyStateRearwards
        "Steady state initialization in rearward direction (chemical potential from fore connector)",
      state
        "Initialization with initial states",
      derivative
        "Initialization with initial derivatives of states") "Choices for initialization of a state."
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(
          coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>
Choices for initialization of a&nbsp;state.
</p>
</html>"));
    type SolutionChoice        = enumeration(
      SolutionPort
        "Chemical solution state from solution port",
      Parameter
        "Chemical solution state from parameter",
      FirstSubstrate
        "Chemical solution state from first substrate")
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(
          coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>
Choices for selection of a&nbsp;chemical solution state.
</p>
</html>"));

    type FirstProductChoice     = enumeration(
      Process
        "First product is defined by process",
      Substance
        "First product definition is a parameter")
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(
        coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>
Choices for selection of a&nbsp;products definition.
</p>
</html>"));

  end Types;

  package Internal "Internal helper functions and models for the undirected themofluid simulation."
    extends Modelica.Icons.InternalPackage;

    function regStep
      "Approximation of a general step, such that the characteristic is continuous and differentiable"
      extends Modelica.Icons.Function;
      input Real x "Abscissa value";
      input Real y1 "Ordinate value for x > 0";
      input Real y2 "Ordinate value for x < 0";
      input Real x_small(min=0) = 1e-5
        "Approximation of step for -x_small <= x <= x_small; x_small >= 0 required";
      output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
    algorithm
      y := smooth(1, if x >  x_small then y1 else
                     if x < -x_small then y2 else
                     if x_small > 0 then (x/x_small)*((x/x_small)^2 - 3)*(y2-y1)/4 + (y1+y2)/2 else (y1+y2)/2);
      annotation(Documentation(revisions="<html>
<ul>
<li><em>April 29, 2008</em>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implemented.</li>
<li><em>August 12, 2008</em>
    by <a href=\"mailto:Michael.Sielemann@dlr.de\">Michael Sielemann</a>:<br>
    Minor modification to cover the limit case <code>x_small -> 0</code> without division by zero.</li>
</ul>
</html>",     info="<html>
<u>
This function is used to approximate the equation
</u>

<blockquote><pre>
y = <strong>if</strong> x &gt; 0 <strong>then</strong> y1 <strong>else</strong> y2;
</pre></blockquote>

<u>
by a smooth characteristic, so that the expression is continuous and differentiable:
</u>

<blockquote><pre>
y = <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> y1 <strong>else</strong>
              <strong>if</strong> x &lt; -x_small <strong>then</strong> y2 <strong>else</strong> f(y1, y2));
</pre></blockquote>

<u>
In the region -x_small &lt; x &lt; x_small a 2nd order polynomial is used
for a smooth transition from y1 to y2.
</u>
</html>"));
    end regStep;

    function regStepState "Apply regStep on State"
      extends Modelica.Icons.Function;

      input Modelica.Units.SI.MolarFlowRate n_flow;
      input Chemical.Interfaces.SubstanceState state_forwards;
      input Chemical.Interfaces.SubstanceState state_rearwards;
      input Modelica.Units.SI.MolarFlowRate n_flow_reg;

      output Chemical.Interfaces.SubstanceState state;

    protected
      Modelica.Units.SI.ChemicalPotential u;
      Modelica.Units.SI.MolarEnthalpy h;

    algorithm
      u := regStep(n_flow, state_forwards.u, state_rearwards.u, n_flow_reg);
      h := regStep(n_flow, state_forwards.h, state_rearwards.h, n_flow_reg);

      state := Chemical.Interfaces.SubstanceState(u=u, h=h);

      annotation (Documentation(info="<html>
<u>This function applies the regStep function to u,T and Xi of a state and creates and returns the resulting state.</u>
</html>"));
    end regStepState;
    annotation (Documentation(info="<html>
<u>Internal helper functions and models for the undirected themofluid simulation.</u>
</html>"));
  end Internal;
end Utilities;
