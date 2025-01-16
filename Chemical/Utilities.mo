within Chemical;
package Utilities
  package Units
    type Inertance =  Real (quantity = "URTInertance",unit = "s2/mol", start=1e-5, nominal=1e-5, min=0) "Inertance of electro-chemical process";
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
  end Types;
end Utilities;
