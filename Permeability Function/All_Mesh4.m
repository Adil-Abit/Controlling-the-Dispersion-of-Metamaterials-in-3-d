function  [Mesh,All_Tets, Vertex]= All_Mesh4()
MY_MESH = load_gmsh2('con0.3(0.12)'); %For a different geometry, choose a different geo file here!
Vertex = MY_MESH.POS;
All_Tets = MY_MESH.TETS;
All_Elements = zeros(length(All_Tets),4);
All_Elements = All_Tets(:,1:4);
Omega = Domain('tetrahedron');
Mesh = MeshTetrahedron(All_Elements, Vertex, 'Omega');
end