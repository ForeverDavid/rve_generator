function ellipsoidIgesOut(x, y, z, outName)
kurve = nrbcirc(0.5,[0,0,0], 0, pi);
oberflaeche_Kugel = nrbrevolve(kurve, [0,0,0], [1,0,0], 2*pi);
oberflaeche_Ellipsoid = nrbtform(oberflaeche_Kugel, vecscale([x, y, z]));
igesout(oberflaeche_Ellipsoid, outName)
display('fertig')