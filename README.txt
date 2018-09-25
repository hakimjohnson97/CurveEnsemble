Curve Ensemble (http://sourceforge.net/projects/curve-ensemble)

Curve Ensemble is a free C++ open-source project, available under the GNU Public license, for constructing, editing, and painting planar curves.  The primary focus is on minimal energy curves, and our implimentation includes (Restricted) Elastic Splines as well as several methods related to parametric cubic splines.

For each minimal energy curve method the following features are available:
   * Robust algorithm to find minimal energy interpolating curve.
   * Flexible node specifications (regular, clamped, free, free-clamped, clamped-free).
   * Multiple curves on the same canvas.
   * Transformation tools (translate, rotate, scale).
   * Elaborate painting capabilities.

There are also other useful methods implemented, like Lines and Circles.

NEW!: Curve Ensemble now allows data balls and can measure tortuosity.

Usage is now fully documented in a 27-page manual (Curve_Ensemble_manual.pdf) which includes three easy-to-follow tutorials as well as very detailed information for the user and the academic.

If you have a curve method you would like to see implemented, please contact
Michael Johnson <yohnson1963@hotmail.com>

****************************************************************************************************************

INSTALLATION

For Windows, the binaries are available at http://sourceforge.net/projects/curve-ensemble
But if you wish to complile it yourself, first download the archive Curve-Ensemble-source-1.2.tar.gz
from http://sourceforge.net/projects/curve-ensemble and then extract it (eg use WinRAR in Windows or the command
       tar -xzf Curve-Ensemble-source-1.2.tar.gz       
in unix).  This will create the directory Curve-Ensemble-source-1.2 containing the source code and other informative files.  To install, follow the directions in the file INSTALL.txt

***************************************************************************************************************

HOW TO USE

Basic instructions:
The basic idea is to specify nodes on the canvas (right mouse-click) and Curve Ensemble will plot a curve through these points using the selected curve method.

------------------------------------------------------------------------------------------------------------
Creating Nodes

To create a node, right click the desired location on the canvas. Individual nodes can be selected by left-clicking on them (you will see a black outline around the selected node). Pressing the delete key will delete the selected node. Note that when you create a node with right-click, it will create a node after the currently selected node. So you can insert a node between two nodes by selecting the previous node and then right clicking. A node can be moved by left-clickling it, holding the left mouse button and dragging it to the desired location.

------------------------------------------------------------------------------------------------------------
Optimization

For Minimal Energy Methods (methods above the line):
Although a curve is drawn after nodes are plotted, it is not optimal. It takes some computing to minimize the energy to make it look smooth. You can presse the optimize button below the canvas to optimize the curve. Fine Optimize will do an even more accurate optimization but is not necessary until you are finalizing the project.

NEW Version 1.2! - Course optimization is now done on the fly. However, fine optimization is still done by clicking the 'fine optimize' button.

------------------------------------------------------------------------------------------------------------
Basic Info

Before we move on to the Curve Editor, it will help to know a bit about how these curves work. It may seem like you specify the points and then a curve is drawn through them and there is nothing else that can be done. But each node also has a direction attached to it - the angle in which the curve comes into the node and out. This direction can either be chosen by the computer (to minimise energy) or it can be chosen by you (where you clamp the direction). So far, the computer has done all the choosing and this is because all the nodes are marked 'Regular'. This makes a very smooth curve with a low bending energy. All this can be changed in the Curve Editor. However, sometimes you may wish to specifiy your own direction at a node to control the shape of the curve. You can do this by changing a node to  'Irregular', clamping the directions and specifying your own direction as a bearing anticlockwise from East. Note there are two directions, the incoming and outgoing. Yes, this means they can be different so the node will have a corner (the angle between the incoming and outgoing directions). 
So far there seems to be two states a node can be, Regular or Irregular with directions clamped. However, when a node is changed to Irregular, it is actually put in the FREE state. This means the node does not actually have a direction and a different shape is made. You can see how this looks by unchecking the Regular box of one of the nodes, the curve may just be a straight line. You may notice that the end nodes (first and last) are actually irregular FREE nodes. This is because it is not a closed curve by defualt. A closed curve can have all the nodes regular. So there are three different states a node can be - REGULAR, IRREGULAR-CLAMPED, IRREGULAR-FREE. Note the incoming and outgoing directions are different so they can have different states.  Read the Curve Editor section to learn how to change the states.

------------------------------------------------------------------------------------------------------------
Curve Editor

This is a separate window that you might have noticed automatically comes up at start-up. It contains all the necessary variables to edit the curve. It is split into two sections: Curve-Specific and Node-Specifc.

Curve-Specific - These affect the curve as a whole

* Closed: This controls whether the curve will form a full loop (join the last node back to the first).

* Curve Method: These are the different basic curve methods that are available, by default Restriced Elastic Spline is used.

Node-Specific - These affect individual nodes (Read Basic Info to understand the states)

* Regularity: Controls if the node is regular or not

* Directions: IF NODE IS IRREGULAR ONLY: Clamp button can be pressed to clamp the direction. The direction can either be typed into the box or the arrow can be edited on the canvas (the arrow will appear when the the direction is clamped). Simply click and hold the left mouse button on the tip of the arrow and drag it around to change the direction. An incoming direction will be displayed as an arrow with a tail coming behind the node, and an outgoing direction will be a conventional arrow going out of the node.

* Corner: This is the angle between the incoming and outgoing directions. This can be changed on a regular node, where the corner will always be maintained. On an irregular node, this is only available if both the incoming and outgoing directions are clamped.  It can be edited but will just changed the incoming direction accordingly. 

* Respect Corner: This is a checkbox next to the Corner box. It is only available when both the incoming and outgoing directions are clamped.  It basically controls whether the corner is held constant when the incoming or outgoing directions are changed. If you move the arrow when respect corner is unchecked, the corner will change. But if it is checked, when one direction changes, the other direction will change accordingly to keep the corner angle constant.  Respect Corner is particularly useful if you want the incoming direction to equal the outgoing direction, but you would like to specify this common direction.  

* Position X and Y: These just provide a more precise way to specify your node positions. You can put two nodes in line with each other if you make the Y-coordinates the same, for example. Note that the coordinates of your mouse cursor are displayed on the bottom left corner of the program, so it is also possible do this with just your mouse.

NOTE: New values in text boxes are only stored after you press [Enter] or navigate to a different field.
WARNING: In text boxes, use [Backspace] rather than [Delete] since the latter will delete the currently selected node.
You should now be able to do everything concerning the technical side of things.

------------------------------------------------------------------------------------------------------------
Multiple Curves

More than one curve can be created on the canvas with all curves being shown simultaneously. This is done using the tab system with each tab representing a different curve. You can click the '+' sign in the top right (just above the canvas) to create a new curve. Curves can also be deleted by clicking the X on any of the tabs. When you create and edit curves, you are only editing one curve at a time which is the tab you have currently selected. To make the curve that is currectly selected stand out, it will be coloured whereas all the other curves will be black. However this can be changed in the Paint Editor.

------------------------------------------------------------------------------------------------------------
Paint Editor

This controls all painting capablities of Curve Ensemble. It is not opened by default on startup and can be accessed by going to View->Paint Editor.  The Curve Editor and Paint Editor windows can be opened and closed via the View menu. It is important to know what the box with 'Use Recommended Display' or 'Use User-Specified Display' means. It basically switches from using paint options that aid with producing the curve (an arrow to control the direction) and paint options that the user specifies. So with recommended display, the selected curve will have blue nodes and red pieces (individual curve pieces corresponding to a node). Any nodes with irregular states and clamped directions will have arrows which can be altered by the mouse. All other curves will be black. This means if you change any colours or Node Style under Recommended Display, you won't see the changes until you switch to User-Specified Display. It is expected that Recommeneded Display is used for the creation and editing of the curves while User-Specified Display is used for the colouring aspects.

These are its features:

* Apply To All Nodes: This is a simple button to quickly apply everything to all nodes (and pieces) on the curve, not just the current node.

* Background Colour: Controls colour of the background.

The rest of the fields will control individual nodes (or pieces). Note that when I refer to a node, I may also be referring to the piece that comes after it just so I won't have to keep saying nodes (or pieces). And if you're wondering what a piece is, it's the part of the curve connecting two nodes; It's a 'piece' of the curve. The rest of the fields will affect the currentely selected node. Most of the fields are self-explanatory but perhaps the last one 'Point Style' needs some explaining. It basically controls what the node looks like, whether it has an arrow showing outgoing direction, or a tail arrow showing incoming direction or maybe you don't want it to be shown at all. Remember, you must change the Recommended Display to User-Specified to see the effects in 'Paint Mode'.

------------------------------------------------------------------------------------------------------------
Zooming and Transformations

Just above the tabs, you may see three things: Normal, Transform and Zoom. Normal will looked pressed in because that is the currently selected one. These are pretty much the different modes the canvas can be in. Pressing any one of them will change the mode. Just remember to press the Normal button when you are done.

* Zoom: This allows you to zoom in and out on the canvas. Left-Cliking will zoom you in, and right clicking will zoom you out. Note it will zoom about the point on the canvas where you clicked.  A double right-click will return you to the original state.

* Transformations:  As soon as this is pressed, a transformation window comes up. This is not a pop-up window, it is to be used with the canvas simultaneously.  Basically you choose your desired transformation with the box at the top and then you have two options: Changing the fields on the window then clicking 'Transform' or Moving a node on the canvas which will carry out the transformation. For example, with 'Translate' you can either translate by a chosen vector (x units the right and y units up) or you can move any one of the nodes on the canvas which will then move the curve as a whole. With more complicated transformations, such as rotate, a green dot will appear on the canvas which will be the pivot of the rotation. After moving the green dot to the desired location, if any of the nodes are moved, it will move such that the whole curve is rotated. For reflections, a line with two dots on the end points can be moved to specify the line of
reflection.

------------------------------------------------------------------------------------------------------------
Saving, Exporting and Importing

The whole project can be saved via File->Save and it can be loaded by File->Open.  The extension .cvp is recommended for projects.

An individual curve can be exported as an ASCII text file (basically a normal .txt file) via File->Export->Text which stores the skeleton of the curve: the nodes, their states and their corresponding directions. It does not save any paint options or whether it is closed. Note this saves only 1 curve, saving as a project will save all the curves. There is also another option of saving to the clipboard which is essentially saving it as a text file but to a file called clipboard.txt. This is usefule for copying curves. For example, let's say you wanted to copy a curve and perhaps do some transformations to it. So on the first curve, you export to clipboard, then you create a new curve (pressing the + sign) and do File->Import->from Clipboard. Now you can do your transformations to one of them.

Another use of Import is if you have created an ASCII text file where each row contains two numbers (the x and y coordinates of a node), then this file can be imported via File->Export->Text.  Curve Ensemble will ignore any preamble at the top of the file, so such files can be easy created in Octave (or Matlab) simply by saving an n x 2 matrix.

The last way of exporting is as an image, you can do this by File->Export->as Image then you choose your desired file type. Curve Ensemble will then ask you to specify your bounding box (you can move the corners of the box on the canvas). Also you can magnify (or increase the resolution of) your image with the scale factor box.  Michael recommends exporting PNG images with scale factor = 2 and then labeling them using GIMP.







