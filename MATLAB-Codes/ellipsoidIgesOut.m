function ellipsoidIgesOut(x, y, z, outName)
kurve = nrbcirc(x/4,[0,0,0], 0, pi);
oberflaeche_Kugel = nrbrevolve(kurve, [0,0,0], [1,0,0], 2*pi);
oberflaeche_Ellipsoid = nrbtform(oberflaeche_Kugel, vecscale([1, y/x, z/x]));
igesout(oberflaeche_Ellipsoid, outName)