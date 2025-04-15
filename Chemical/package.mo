within ;
package Chemical "Chemical pathways"

 extends Modelica.Icons.Package;

  annotation (
preferredView="info",
version="2.0.0-alpha",
versionDate="2025-04-15",
dateModified = "2025-04-15 15:45:41Z",
conversion(
  from(version="1.4.1", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.4_to_2.0.mos",
        to="2.0.0"),
  from(version="1.4.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.4_to_2.0.mos",
        to="2.0.0"),
  from(version="1.3.1", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_2.0.mos",
        to="2.0.0"),
  from(version="1.3.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_2.0.mos",
        to="2.0.0"),
  from(version="1.2.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.3_to_2.0.mos",
        to="2.0.0"),
  from(version="1.1.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.1_to_2.0.mos",
        to="2.0.0"),
  from(version="1.0.0", script="modelica://Chemical/Resources/Scripts/Dymola/ConvertChemical_from_1.0_to_2.0.mos",
        to="2.0.0")),
      uses( Modelica(version="4.0.0")),
  Documentation(revisions="<html>
<p>Copyright (c) 2025, Marek Matej&aacute;k, Ph.D. </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info="<html>
<p>During each electro-chemical process an <a href=\"modelica://Chemical.Obsolete.Components.Substance\">electro-chemical potential</a> of the substances is equilibrating and all thermodynamical properties of the homogenous chemical solutions are evaluated. </p>
<p>Processes: chemical reactions, gas dissolution, diffusion, membrane transports, osmotic fluxes, electrochemical cells, electrodes, ..</p>
<p>Please see the <a href=\"modelica://Chemical.UsersGuide.Overview\">overview</a>.</p>
</html>"));
end Chemical;
