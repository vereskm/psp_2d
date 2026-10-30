function [shapes] = getShapes()
    %GETSHAPES RETURNS THE AVAILABLE SHAPES. EACH SHOULD HAVE AN INSIDE


    %Edit in the same way to get more shapes
    shapes.rod= @geometry.shapes.insideRod;
    shapes.waveg=@geometry.shapes.insideWaveguide;
    shapes.ring=@geometry.shapes.insideRing;


end