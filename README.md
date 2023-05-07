Download Link: https://assignmentchef.com/product/solved-fem-project-division-of-solid-mechanics
<br>



The task is to write a finite element program and use the program to analyze temperature and stress distributions in an optical system. The problem should be solved with Matlab along with suitable subroutines included in the CALFEM toolbox. A well structured report that presents the findings should be written.

<h1>Problem description</h1>

On February 18, 2021, the NASA rover named Perseverance landed on Mars with a mission to search for signs of ancient life and collect material samples. To perform the scientific observations, and to navigate the treacherous landscape of Mars, it has been equipped with no less than 23 cameras. Obviously, the functionality of these cameras are of uttermost importance for the mission to end successfully.

Due to the thin atmosphere, the temperatures are known to vary considerably during a Martian day. Unfortunately, the temperature differences will develop thermal stresses in the lenses inside the cameras, as these are constructed from multiple types of materials and might be geometrically constrained. The associated deformation might interfere with the optical properties.

This must be taken into account when designing a lens on Perseverance. To this end, a finite element analysis is performed. Our simplified model of a lens consists of two different materials. The casing is made from a titanium (Ti) alloy and the lenses are made of glass (Gl). A cross section of the lens that should be analyzed is depicted in Fig. 2.

Figure 1: Sketch of a cross section of the lens. All dimensions are in cm.

Since the two lenses are curved, these must be meshed using two ellipses, both with semi-axes 0<em>.</em>1 cm in the <em>x</em>-direction and 0<em>.</em>25 cm in the <em>y</em>-direction, cf. Fig. 2.

Figure 2: Drawing of the lens in PDE tool.

The material data are provided in Table 1.

<table width="392">

 <tbody>

  <tr>

   <td width="225"> </td>

   <td width="116">Titanitum alloy</td>

   <td width="52">Glass</td>

  </tr>

  <tr>

   <td width="225">Young’s modulus, <em>E </em>[GPa]</td>

   <td width="116">110</td>

   <td width="52">67</td>

  </tr>

  <tr>

   <td width="225">Poisson’s ratio, <em>ν </em>[-]</td>

   <td width="116">0.34</td>

   <td width="52">0.2</td>

  </tr>

  <tr>

   <td width="225">Expansion coefficient, <em>α </em>[1<em>/K</em>]</td>

   <td width="116">9<em>.</em>4 · 10<sup>−6</sup></td>

   <td width="52">7 · 10<sup>−6</sup></td>

  </tr>

  <tr>

   <td width="225">Density, <em>ρ </em>[kg/m<sup>3</sup>]</td>

   <td width="116">4620</td>

   <td width="52">3860</td>

  </tr>

  <tr>

   <td width="225">Specific heat, <em>c<sub>p </sub></em>[J/(kg K)]</td>

   <td width="116">523</td>

   <td width="52">670</td>

  </tr>

  <tr>

   <td width="225">Thermal conductivity, <em>k </em>[W/(m K)]</td>

   <td width="116">17</td>

   <td width="52">0.8</td>

  </tr>

 </tbody>

</table>

Table 1: Material data

The lens protrudes 1 cm out-of-plane (i.e. thickness). It is mounted in the camera such that the top and bottom boundaries are fixated in the <em>y</em>-direction, i.e. <em>u<sub>y </sub></em>= 0, whereas the rightmost boundary is constrained in the <em>x</em>-direction, i.e. <em>u<sub>x </sub></em>= 0. Fig. 2 illustrates the boundary conditions. Boundaries without displacement boundary conditions are traction free. Plane strain conditions are assumed to hold. Newton convection is present along parts of the boundary of the lens, cf. Fig. 2, i.e., <em>q<sub>n </sub></em>= <em>α<sub>c</sub></em>(<em>T </em>− <em>T</em>∞) on the leftmost side and <em>q<sub>n </sub></em>= <em>α<sub>c</sub></em>(<em>T </em>− <em>T<sub>c</sub></em>) on the rightmost side, where <em>α<sub>c </sub></em>= 100 W/(m<sup>2 </sup>K). The volumes in between the lenses are assumed to contain vacuum, making the internal boundaries thermally insulated. The remaining external boundaries are also thermally insulated. The structure is stress free at <em>T</em>0 = 20 <sup>◦</sup>C. If possible you should use symmetry conditions.

<h2>Problem formulation</h2>

<ol>

 <li>Find the stationary temperature distribution for daytime conditions, i.e. when the lens isexposed to ambient air with temperature of <em>T</em>∞ = 40 <sup>◦</sup>C as well as a temperature control system that ensures <em>T<sub>c </sub></em>= 20 <sup>◦</sup>C along the rightmost boundary. Do the same for night time conditions, assuming an ambient temperature of <em>T</em>∞ = −96 <sup>◦</sup>C and again <em>T<sub>c </sub></em>= 20 <sup>◦</sup> Compute the maximum temperature for both stationary conditions.</li>

 <li>Determine the transient temperature evolution during daytime, assuming that the initialconditions coincide with the stationary night conditions. Also, compute the temperature evolution during night assuming that the initial conditions are given by the stationary day conditions. Present the temperature distributions at meaningful time instants.</li>

 <li>Solve the mechanical problem. Determine and present the effective von Mises stress fieldfor both night and day, i.e. using the same conditions as in task a). Where are high stress concentrations found?</li>

 <li>In order for the lenses to maintain their intended optical properties their curvature should remain as close as possible to their original design during temperature shifts. One way to quantify the displacement is to compute the square sum of the displacements in the lens. Determine this value for the leftmost lens subject to the thermal conditions in task a) and present also the displacement field (deformation pattern) when stationary temperature is reached during day/night.</li>

</ol>

Hint 1: The von Mises stress is defined as:

<em>σ</em><em>eff </em>= q<em>σ</em><em>xx</em>2 + <em>σ</em><em>yy</em>2 + <em>σ</em><em>zz</em>2 − <em>σ</em><em>xx</em><em>σ</em><em>yy </em>− <em>σ</em><em>xx</em><em>σ</em><em>zz </em>− <em>σ</em><em>yy</em><em>σ</em><em>zz </em>+ 3<em>τ</em><em>xy</em>2 + 3<em>τ</em><em>xz</em>2 + 3<em>τ</em><em>yz</em>2 Hint 2: The displacement magnitude squared for the leftmost lens Ω<em><sub>lens </sub></em>can be found as

u<em><sup>T</sup></em>u<em>dV </em>= a<em><sup>T </sup></em>Z                   N<em><sup>T</sup></em>N<em>dV </em>a = a<em><sup>T</sup></em>Ta

Ω<em>lens                                                  </em>Ω<em>lens</em>

<h1>Procedure</h1>

A fully implicit time integration scheme should be used. Note that the element function for forming, C<em><sup>e</sup></em>, is available in Canvas, and that T<em><sup>e </sup></em>is obtained by modifying the routine that computes C<em><sup>e</sup></em>. A suitable element is the linear triangular element. As a start point, you have the strong formulations of the heat and the mechanical problems.

The contour plots of the stress distribution are based on the stress at the nodal points. The stress of the nodal points can be approximated by taking the mean value of the stresses in the elements connected to a node. The following Matlab code can be useful:

for i=1:size(coord,1)

[c0,c1]=find(edof(:,2:4)==i);

Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1); end

where Seff_nod and Seff_el is the von Mises effective stress at the nodal points and in the elements, respectively. Note that here the topology matrix (edof) is associated with the temperature problem.

<h1>Report</h1>

A fundamental ingredient in all research is that it should be possible to regenerate the results obtained based on the report. In the present situation this implies that the appended matlab code should only be considered as supporting material. Moreover, note that one variable for grading the report is the structure of the computer code, i.e. you should choose suitable names for variables. A suitable structure for the report is:

<ul>

 <li>Introduction: Description of the problem, geometry and boundary conditions.</li>

 <li>Procedure: How the problems are solved (weak formulation, application of boundary conditions, treatment of thermal strains and so on). Derivation of the FE formulation and other important theoretical aspects should be included here. Note that you are encouraged to make references to textbooks. Though, it is important to carefully present all calculations that are not available in the literature.</li>

 <li>Results: Present the results in illustrative figures and/or tables. Note that the results should be commented such that the reader can not misunderstand the results (correct labels, units, figure texts etc.)</li>

 <li>Discussion: A discussion of the results, their meaning and significance. You might want to discuss sources of errors and accuracy in this section.</li>

 <li>Computer Code: Note that the code should be easy to follow and all declared variables should have intuitive names and be well commented.</li>

</ul>

A well structured should be handed to the Division of Solid Mechanics no later than 2021-05-26 16.00. Matlab/CALFEM files (appendix) should be well commented. The reader of the report is assumed to have the same knowledge level as the author. If the report contains major programming or theoretical errors, the report is returned in order to be corrected. It is possible to obtain up to 5 points which are added to the points obtained at the exam in June 2021. The assignment should be approved no later than 2021-06-09. In addition to your report in PDF format, you should also attach your m-files in the email. Note that the bonus points obtained is only valid for the examination in June 2021.

<h1>1           Generating mesh</h1>

In this assignment you should create your own mesh using the built-in PDEtool in MATLAB. This tutorial considers the example geometry seen in figure 3, which of course must be changed for your specific application.

PDEtool can be used directly to create simple geometries with little effort. To start just type pdetool into the MATLAB terminal to open PDEtool user interface. Although the tool is originally designed for solving partial differential equations only the meshing aspect are of interest in this assignment.

Figure 3: Dog toy geometry.

A step by step tutorial of how to use some of the basic features of PDEtool are provided here. Before starting to draw it is often convenient to activate grid and change the axis scales. This is done using &lt;Options&gt; menu. To build our geometry simple geometries are added one by one. The basic shapes used in PDEtool are rectangles, ellipses and polygons. Basic drawing tools can be found on the toolbar (see figure 4, note polygons are not necessary to create the geometry for the lab).

Figure 4: Basic tools for drawing in PDEtool.

Draw a rectangle: Create a rectangle by using the rectangle tool (see figure 4) and simply hold left mouse button and move the mouse cursor to obtain the desired size. To adjust the size and position double click on the rectangle and a dialog will appear where coordinates and dimensions can be provided, see figure 5.

Figure 5: Rectangle object dialog.

Draw an ellipse: Create a ellipse by using the draw ellipse tool (see figure 4) and simply hold left mouse button and move the mouse cursor to obtain the desired size. Double click to adjust position and size in similar fashion as for rectangles.

Determining combined geometry: In the box below shape-buttons the present features are presented by name (which may be edited using object dialog). The plus sign indicates that geometries are added and minus signs indicate that they are cut out. Because of this the order is important i.e R1+E1-E2 =6 R1-E2+E1 in general. Below figues showing drawn features and the resulting body after combining by grouping and changing signs in the set formula row (see figure 6, 8 and 7). Note that to see the resulting geometry we need to go to either meshing or boundary view (see toolbar).

Figure 6: Set of combined shapes.

When the necessary features are added go to &lt;Boundary&gt; / &lt;Boundary mode&gt; to see which are the boundaries of your geometry. The red arrows define the main borders and the grey contours represent the so called subdomains. For instance, subdomains appear where features overlap. To get a uniform mesh it is convenient to remove the subdomains. To do so use &lt;Boundary&gt; / &lt;Remove all Subdomain Borders&gt;. If the mesh consist of different natural domains such as different materials it is favourable to keep the subdomains.

Figure 7: Resulting geometry seen in boundary mode (ctrl+B) after all subdomains have been removed.

To obtain the geometry provided in figure 7 the formula used is given in figure 8 and the object names is seen in figure 6.

Figure 8: Formula determining active areas.

To generate a mesh simply press mesh tool when the geometry seen in boundary mode is satisfactory (see figure 4). To refine the mesh, i.e. create more elements, press refine mesh until enough elements are generated (remember that finer mesh means higher computational cost but better accuracy). Note that a high number of elements increases the computational time.

Figure 9: Mesh after two refinements.

When you are content with your mesh use &lt;Mesh&gt; / &lt;Export mesh&gt; to save the topology matrices associated with the current mesh (note you have to save the work separately to keep geometries etc. since &lt;Export mesh&gt; will only save matrices). The topology matrices are p,e,t (points, edges, triangles). The exported matrices will directly appear in the active MATLAB workspace. It is strongly advised to save them directly in a .mat file (mark variables in workspace, right click and save) so it can be loaded in the scrips you write. From p,e,t the CALFEM quantities edof, dof, coord etc can be extracted. On the course web page for FEM FAQ are some instructions of how this is done.

<h1>2           General tips</h1>

The output data from PDEtool is in the form p,e,t where p is points, it is a 2 row matrix where the first row is the x-coordinates and the second row is the y-coordinates. There are one column for each node i.e. in CALFEM notation coord = p’.

The output t is triangles and it has dimensions 4×nelm. The three first rows are the node numbers associated with each element (three nodes for linear triangular elements) while the fourth row is its subdomain. The subdomain is great to use when you want to associate specific material properties to each element. Rows 1,2 and 3 can be used to calculate our edof matrices, as an example

<table width="454">

 <tbody>

  <tr>

   <td width="265">enod=t(1:3,:)’;</td>

   <td width="189">% nodes of elements</td>

  </tr>

  <tr>

   <td width="265">nelm=size(enod,1);</td>

   <td width="189">% number of elements</td>

  </tr>

  <tr>

   <td width="265">nnod=size(coord,1);</td>

   <td width="189">% number of nodes</td>

  </tr>

  <tr>

   <td width="265">dof=(1:nnod)’;</td>

   <td width="189">% dof number is node number</td>

  </tr>

 </tbody>

</table>

dof_S=[(1:nnod)’,(nnod+1:2*nnod)’]; % give each dof a number for ie=1:nelm edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)]; edof(ie,:)=[ie,enod(ie,:)];

end

The output e is edges and contains among other things node-pairs belonging to a certain boundary segment. Thus if you activate the option in PDEtool, &lt;Show Edge Labels&gt; and &lt;Show Subdomain Labels&gt; in boundary mode you can see how the different segments are numbered in e and subdomains in t. For our intents and purposes the rows 1,2 and 5 are the relevant ones. Row 1 and 2 includes the node numbers of the element-segment and row 5 is the edge label, i.e. which boundary segment it belongs to. For instance a list of all the convective boundaries for the heat problems can be found from

% Check which segments that should have convections er = e([1 2 5],:); % Reduced e

%conv_segments = [10 11 12]; % Choosen boundary segments edges_conv = []; for i = 1:size(er,2) if ismember(er(3,i),conv_segments) edges_conv = [edges_conv er(1:2,i)];

end

end

edges_conv now contains the node numbers of the node pairs that belongs to the convective part of the global boundary. From this the length of each segment <em>L </em>can be determined when calculating the contributions to boundary force vector and the stiffness matrix. Figure 10 shows an example of how the segment numbering and subdomains can look like.

Figure 10: Lens geometry with subdomains (red numbers) and edge labels (black numbers). Note that the meshing will be constrained by the edges (grey lines) so that an element is wholly within a single subdomain. Boundary segment 1 and 14 in this case corresponds to where the boundary conditions are changing at the outer boundary, e.g. 14 is isolated whereas 1 has convection.

Having made a mesh in PDE it is trivial to refine the mesh so it is advantageous to use a small mesh to test your code since it is much faster. When you think it is working use a finer mesh to get better resolved results.

From experience the CALFEM function assem and insert are very slow. The following code rows does the same thing:

% Kt = assem(edof(el,:),Kt,Kte); indx = edof(el,2:end);

Kt(indx,indx) = Kt(indx,indx)+Kte;

% f = insert(edof(el,:),f,fe); indx = edof(el,2:end); f(indx) = f(indx) + fe;

It is always nice to present nice plots, here are some short tips for plotting.

For instance the Matlab function patch can be used to generate field plots, e.g. patch(ex’,ey’,eT’) where eT is the element temperatures obtained from eT=extract(edof,T) with T as the nodal temperatures. To plot the whole component (both sides of the symmetry cut) you can therefore write:

patch(ex’,ey’,eT’) hold on patch(-ex’,ey’,eT’)

If you use a fine mesh, the mesh lines can sometimes obscure the actual result. To remove the mesh lines from the plotting the option ’EdgeColor’ can be turned off (note that somewhere in the report the mesh used for the results should be shown), i.e. patch(ex’,ey’,eT’,’EdgeColor’,’none’).

Don’t forget to turn on the colorbar and choose an illustrative color scale, e.g. colormap(hot) and set the axis scales to the correct proportions. Example:

figure() patch(ex’,ey’,eT’,’EdgeColor’,’none’) title(’Temperature distribution [C]’) colormap(hot); colorbar; xlabel(’x-position [m]’) ylabel(’y-position [m]’) axis equal

In order to compare several plots with each other it is recommended to use the same scale. By setting option caxis([Tmin Tmax]) where Tmin and Tmax can be calculated or chosen.

To compare the deformation patterns it is preferred they are plotted on top of each other. This can be done in several ways but if patch is used the opacity can be set by the option ’FaceAlpha’ where 1 is completely opaque and 0 completely transparent. Example:

% Calculate displaced coordinates

mag = 100; % Magnification (due to small deformations) exd = ex + mag*ed(:,1:2:end); eyd = ey + mag*ed(:,2:2:end);

figure()

patch(ex’,ey’,[0 0 0],’EdgeColor’,’none’,’FaceAlpha’,0.3) hold on

patch(exd’,eyd’,[0 0 0],’FaceAlpha’,0.3) axis equal

title(’Displacement field [Magnitude enhancement 100]’)

where the third argument is a color triplet, i.e. [0 0 0] is black, so with ’FaceAlpha’=0.3 the result will be a grey-scale plot. You might have to change the magnitude to get nice  plots