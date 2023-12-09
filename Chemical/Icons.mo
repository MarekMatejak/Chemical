within Chemical;
package Icons "Icons for chemical models"
  //extends Modelica.Icons.IconsPackage;
  extends Modelica.Icons.Package;

  partial class Diffusion

    annotation (Icon(graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://Chemical/Resources/Icons/diffusion.png")}));

  end Diffusion;

  class Substance

      annotation ( Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
          graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://Chemical/Resources/Icons/Substance.png")}));
  end Substance;

  class Speciation

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
          graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://Chemical/Resources/Icons/Speciation.png")}));
  end Speciation;

  class GasSolubility

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Bitmap(extent={{-100,-100},{100,
              100}},
              fileName="modelica://Chemical/Resources/Icons/GasSolubility.png")}));
  end GasSolubility;

  class Membrane

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Bitmap(extent={{-100,100},{100,-100}},
              fileName="modelica://Chemical/Resources/Icons/membrane.png")}));
  end Membrane;

  class EnzymeKinetics

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Bitmap(extent={{-80,-26},{86,84}},
              fileName="modelica://Chemical/Resources/Icons/EnzymeKinetics.png")}));
  end EnzymeKinetics;

  class Solution

      annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={
        Line(
          points={{-98,90},{-94,96},{-84,98},{84,98},{96,96},{100,92},{98,86},
              {94,80},{94,80},{94,-92},{94,-92},{94,-96},{92,-100},{88,-100},
              {84,-100},{-84,-100},{-88,-100},{-92,-100},{-94,-96},{-94,-92},
              {-94,24},{-94,78},{-94,80},{-98,90}},
          color={127,0,127},
          smooth=Smooth.Bezier,
          pattern=LinePattern.Dot,
          thickness=0.5)}));
  end Solution;

  class Buffer

      annotation ( Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
          graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://Chemical/Resources/Icons/buffer.png")}));
  end Buffer;

  class ElectronTransfer

      annotation ( Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
          graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://Chemical/Resources/Icons/electron.png")}));
  end ElectronTransfer;
  annotation (Documentation(revisions=""));
end Icons;
