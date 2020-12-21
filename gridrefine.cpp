#include "header.h"

std::vector<std::string> tagNames;

void Refine(Mesh *m)
{
    double cnt[3];
    ElementArray<Node> edge_nodes(m,2); // auxiliary node pair to introduce new edge
    ElementArray<Edge> face_edges(m,4); // my
    MarkerType node_cedge = m->CreateMarker(); // nodes in edge centers
    MarkerType new_edge   = m->CreateMarker(); // my
    MarkerType node_cface = m->CreateMarker(); // nodes in face centers

    m->BeginModification();

    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        if(icell->New())
            continue;

        ElementArray<Face> face(m);

        ///////////////////////////////////////////////////////
        /// Split edges of the cell
        ElementArray<Edge> cellEdges = icell->getEdges();
        for(auto iedge = cellEdges.begin(); iedge != cellEdges.end(); iedge++){
            if(iedge->New())
                continue;
            for(int k = 0; k < 3; ++k)
                cnt[k] = (iedge->getBeg().Coords()[k]+iedge->getEnd().Coords()[k])*0.5;
            Node n = m->CreateNode(cnt);
            n.SetMarker(node_cedge);
            Edge::SplitEdge(iedge->self(),ElementArray<Node>(m,1,n.GetHandle()),0);
        }


        ///////////////////////////////////////////////////////
        /// Split faces of the cell
        ElementArray<Face> cellFaces = icell->getFaces();
        for(auto iface = cellFaces.begin(); iface != cellFaces.end(); iface++){
            if(iface->New())
                continue;

            ElementArray<Node> new_nodes = iface->getNodes(node_cedge);
            ElementArray<Edge> new_edges(m,new_nodes.size());
            if( new_nodes.size() == 3)
            {
                std::cout << "new_nodes.size() == 3" << std::endl;
                exit(3);
            }
            else
            {
                printf("size = %d\n", new_nodes.size());
                cnt[0] = cnt[1] = cnt[2] = 0;
                for(int k = 0; k < new_nodes.size(); ++k)
                {
                    cnt[0] += new_nodes[k].Coords()[0];
                    cnt[1] += new_nodes[k].Coords()[1];
                    cnt[2] += new_nodes[k].Coords()[2];
                }
                cnt[0] /= new_nodes.size();
                cnt[1] /= new_nodes.size();
                cnt[2] /= new_nodes.size();

                // Create a new node in the face center
                Node n = m->CreateNode(cnt);
                n->SetMarker(node_cface);

                // Create a set of edges adjacent to the new node
                edge_nodes[0] = n;
                for(int k = 0; k < new_nodes.size(); ++k)
                {
                    edge_nodes[1] = new_nodes[k];
                    new_edges[k] = m->CreateEdge(edge_nodes).first;
                    new_edges[k].SetMarker(new_edge);
                }
            }
            auto splitf = Face::SplitFace(iface->self(),new_edges,0);
//            for(auto iff = splitf.begin(); iff != splitf.end(); iff++)
//                face.push_back(iff->getAsFace());
        }

        ElementArray<Node> nodes = icell->getNodes();

        // Determine cell boundaries
        double maxcnt[3], mincnt[3];
        std::fill_n(maxcnt, std::numeric_limits<double>::min(), 3);
        std::fill_n(mincnt, std::numeric_limits<double>::max(), 3);
        for(auto inode = nodes.begin(); inode != nodes.end(); inode++){
            if(inode->GetMarker(node_cedge) || inode->GetMarker(node_cface))
                continue;

            inode->Barycenter(cnt);
            if(cnt[0] > maxcnt[0] && cnt[1] > maxcnt[1] && cnt[2] > maxcnt[2]){
                maxcnt[0] = cnt[0];
                maxcnt[1] = cnt[1];
                maxcnt[2] = cnt[2];
            }
            if(cnt[0] < mincnt[0] && cnt[1] < mincnt[1] && cnt[2] < mincnt[2]){
                mincnt[0] = cnt[0];
                mincnt[1] = cnt[1];
                mincnt[2] = cnt[2];
            }
        }

        double xmin = mincnt[0];
        double ymin = mincnt[1];
        double zmin = mincnt[2];
        double xmax = maxcnt[0];
        double ymax = maxcnt[1];
        double zmax = maxcnt[2];
        double dx = (xmax - xmin)/2.;
        double dy = (ymax - ymin)/2.;
        double dz = (zmax - zmin)/2.;

        //printf("dx = %e, dy = %e, dz = %e\n", dx, dy, dz);

        //printf("Cell %d: [%lf %lf]  x  [%lf %lf]  x  [%lf %lf]\n", icell->LocalID(), xmin, xmax, ymin, ymax, zmin, zmax);

        Node nnodes[27];
        // Determine nodes 0, 2, 6, 8, 18, 20, 24, 26
        double eps = 1e-7;
        for(auto inode = nodes.begin(); inode != nodes.end(); inode++){
            inode->Barycenter(cnt);
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[0] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[2] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[6] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[8] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[18] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[20] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[24] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[26] = inode->getAsNode();
            }
        }

        // Manual specification for the rest
        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin;
        nnodes[1] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[3] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[4] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[5] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin;
        nnodes[7] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[9] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[10] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[11] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[12] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[13] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[14] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[15] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[16] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[17] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmax;
        nnodes[19] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[21] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[22] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[23] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmax;
        nnodes[25] = m->CreateNode(cnt);


        std::cout << "Created nodes for cell " << icell->LocalID() << std::endl;

        ElementArray<Node> nd(m,4); // nodes to form a face

        nd[0] = nnodes[1];
        nd[1] = nnodes[4];
        nd[2] = nnodes[13];
        nd[3] = nnodes[10];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[4];
        nd[1] = nnodes[7];
        nd[2] = nnodes[16];
        nd[3] = nnodes[13];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[10];
        nd[1] = nnodes[13];
        nd[2] = nnodes[22];
        nd[3] = nnodes[19];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[13];
        nd[1] = nnodes[16];
        nd[2] = nnodes[25];
        nd[3] = nnodes[22];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[9];
        nd[1] = nnodes[10];
        nd[2] = nnodes[13];
        nd[3] = nnodes[12];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[10];
        nd[1] = nnodes[11];
        nd[2] = nnodes[14];
        nd[3] = nnodes[13];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[12];
        nd[1] = nnodes[13];
        nd[2] = nnodes[16];
        nd[3] = nnodes[15];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[13];
        nd[1] = nnodes[14];
        nd[2] = nnodes[16];
        nd[3] = nnodes[17];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[3];
        nd[1] = nnodes[4];
        nd[2] = nnodes[13];
        nd[3] = nnodes[12];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[4];
        nd[1] = nnodes[5];
        nd[2] = nnodes[14];
        nd[3] = nnodes[13];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[12];
        nd[1] = nnodes[13];
        nd[2] = nnodes[22];
        nd[3] = nnodes[21];
        face.push_back(m->CreateFace(nd).first);

        nd[0] = nnodes[13];
        nd[1] = nnodes[14];
        nd[2] = nnodes[23];
        nd[3] = nnodes[22];
        face.push_back(m->CreateFace(nd).first);

        printf("Gonna split with %zu faces\n", face.size());

        // Now introduce 36 faces to split the cell
        /*ElementArray<Face> face(m, 36);
        ElementArray<Node> nd(m,2); // nodes to form an edge
        ElementArray<Edge> ed(m,4); // edges to form a face

        nd[0] = nnodes[0];
        nd[1] = nnodes[1];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[1];
        nd[1] = nnodes[4];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[4];
        nd[1] = nnodes[3];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[3];
        nd[1] = nnodes[0];
        ed.push_back(m->CreateEdge(nd).first);
        face[0] = m->CreateFace(ed).first;

        nd[0] = nnodes[1];
        nd[1] = nnodes[2];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[2];
        nd[1] = nnodes[5];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[5];
        nd[1] = nnodes[4];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[4];
        nd[1] = nnodes[1];
        ed.push_back(m->CreateEdge(nd).first);
        face[1] = m->CreateFace(ed).first;

        nd[0] = nnodes[3];
        nd[1] = nnodes[4];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[4];
        nd[1] = nnodes[7];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[7];
        nd[1] = nnodes[6];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[6];
        nd[1] = nnodes[3];
        ed.push_back(m->CreateEdge(nd).first);
        face[2] = m->CreateFace(ed).first;

        nd[0] = nnodes[4];
        nd[1] = nnodes[5];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[5];
        nd[1] = nnodes[8];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[8];
        nd[1] = nnodes[7];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[7];
        nd[1] = nnodes[5];
        ed.push_back(m->CreateEdge(nd).first);
        face[3] = m->CreateFace(ed).first;

        nd[0] = nnodes[9];
        nd[1] = nnodes[10];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[10];
        nd[1] = nnodes[13];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[13];
        nd[1] = nnodes[12];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[12];
        nd[1] = nnodes[9];
        ed.push_back(m->CreateEdge(nd).first);
        face[4] = m->CreateFace(ed).first;

        nd[0] = nnodes[10];
        nd[1] = nnodes[11];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[11];
        nd[1] = nnodes[14];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[14];
        nd[1] = nnodes[13];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[13];
        nd[1] = nnodes[10];
        ed.push_back(m->CreateEdge(nd).first);
        face[5] = m->CreateFace(ed).first;

        nd[0] = nnodes[12];
        nd[1] = nnodes[13];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[13];
        nd[1] = nnodes[16];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[16];
        nd[1] = nnodes[15];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[15];
        nd[1] = nnodes[12];
        ed.push_back(m->CreateEdge(nd).first);
        face[6] = m->CreateFace(ed).first;

        nd[0] = nnodes[13];
        nd[1] = nnodes[14];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[14];
        nd[1] = nnodes[17];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[17];
        nd[1] = nnodes[16];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[16];
        nd[1] = nnodes[13];
        ed.push_back(m->CreateEdge(nd).first);
        face[7] = m->CreateFace(ed).first;

        nd[0] = nnodes[18];
        nd[1] = nnodes[19];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[19];
        nd[1] = nnodes[22];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[22];
        nd[1] = nnodes[21];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[21];
        nd[1] = nnodes[18];
        ed.push_back(m->CreateEdge(nd).first);
        face[8] = m->CreateFace(ed).first;

        nd[0] = nnodes[19];
        nd[1] = nnodes[20];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[20];
        nd[1] = nnodes[23];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[23];
        nd[1] = nnodes[22];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[22];
        nd[1] = nnodes[19];
        ed.push_back(m->CreateEdge(nd).first);
        face[9] = m->CreateFace(ed).first;

        nd[0] = nnodes[21];
        nd[1] = nnodes[22];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[22];
        nd[1] = nnodes[25];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[25];
        nd[1] = nnodes[24];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[24];
        nd[1] = nnodes[21];
        ed.push_back(m->CreateEdge(nd).first);
        face[10] = m->CreateFace(ed).first;

        nd[0] = nnodes[22];
        nd[1] = nnodes[23];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[23];
        nd[1] = nnodes[26];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[26];
        nd[1] = nnodes[25];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[25];
        nd[1] = nnodes[22];
        ed.push_back(m->CreateEdge(nd).first);
        face[11] = m->CreateFace(ed).first;

        // starting vertical faces
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;

        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed[0] = m->CreateEdge(nd).first;
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        nd[0] = nnodes[];
        nd[1] = nnodes[];
        ed.push_back(m->CreateEdge(nd).first);
        face[] = m->CreateFace(ed).first;*/

        auto split = Cell::SplitCell(icell->self(), face, 0);
        printf("Split: got %zu cells\n", split.size());
    }

    m->EndModification();
    m->ReleaseMarker(node_cedge,NODE);

    m->AssignGlobalID(CELL|FACE);

    m->Save("refined.vtk");
}

void Refine1(Mesh *m)
{
    double cnt[3];
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        if(icell->New())
            continue;
        ElementArray<Node> nodes = icell->getNodes();

        // Determine cell boundaries
        double maxcnt[3], mincnt[3];
        std::fill_n(maxcnt, std::numeric_limits<double>::min(), 3);
        std::fill_n(mincnt, std::numeric_limits<double>::max(), 3);
        for(auto inode = nodes.begin(); inode != nodes.end(); inode++){

            inode->Barycenter(cnt);
            if(cnt[0] > maxcnt[0] && cnt[1] > maxcnt[1] && cnt[2] > maxcnt[2]){
                maxcnt[0] = cnt[0];
                maxcnt[1] = cnt[1];
                maxcnt[2] = cnt[2];
            }
            if(cnt[0] < mincnt[0] && cnt[1] < mincnt[1] && cnt[2] < mincnt[2]){
                mincnt[0] = cnt[0];
                mincnt[1] = cnt[1];
                mincnt[2] = cnt[2];
            }
        }

        double xmin = mincnt[0];
        double ymin = mincnt[1];
        double zmin = mincnt[2];
        double xmax = maxcnt[0];
        double ymax = maxcnt[1];
        double zmax = maxcnt[2];
        double dx = (xmax - xmin)/2.;
        double dy = (ymax - ymin)/2.;
        double dz = (zmax - zmin)/2.;

        printf("dx = %e, dy = %e, dz = %e\n", dx, dy, dz);

        printf("Cell %d: [%lf %lf]  x  [%lf %lf]  x  [%lf %lf]\n", icell->LocalID(), xmin, xmax, ymin, ymax, zmin, zmax);

        Node nnodes[27];
        // Determine nodes 0, 2, 6, 8, 18, 20, 24, 26
        double eps = 1e-7;
        for(auto inode = nodes.begin(); inode != nodes.end(); inode++){
            inode->Barycenter(cnt);
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[0] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[2] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[6] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmin) < eps){
                nnodes[8] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[18] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymin) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[20] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmin) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[24] = inode->getAsNode();
            }
            if(fabs(cnt[0] - xmax) < eps && fabs(cnt[0] - ymax) < eps && fabs(cnt[0] - zmax) < eps){
                nnodes[26] = inode->getAsNode();
            }
        }

        // Manual specification for the rest
        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin;
        nnodes[1] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[3] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[4] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nnodes[5] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin;
        nnodes[7] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[9] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[10] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nnodes[11] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[12] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[13] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nnodes[14] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[15] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[16] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nnodes[17] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmax;
        nnodes[19] = m->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[21] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[22] = m->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nnodes[23] = m->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmax;
        nnodes[25] = m->CreateNode(cnt);

        //std::cout << "Created nodes for cell " << icell->LocalID() << std::endl;

        ElementArray<Face> facess = icell->getFaces();
        //double eps = 1e-8;
        for(auto iface = facess.begin(); iface != facess.end(); iface++){
            iface->Normal(cnt);
            for(int i = 0; i < 3; i++)
                cnt[i] /= iface->Area();

            if(fabs(cnt[0]) < eps && fabs(cnt[1]+1.0) < eps && fabs(cnt[2]) < eps){
                printf("Found face, nor = %e %e %e\n", cnt[0], cnt[1], cnt[2]);

                ElementArray<Edge> splittingEdges(m);
                ElementArray<Node> n(m,2);
                n[0] = nnodes[10];
                n[1] = nnodes[1];
                splittingEdges.push_back(m->CreateEdge(n).first);
                n[0] = nnodes[10];
                n[1] = nnodes[9];
                splittingEdges.push_back(m->CreateEdge(n).first);
                n[0] = nnodes[10];
                n[1] = nnodes[19];
                splittingEdges.push_back(m->CreateEdge(n).first);
                n[0] = nnodes[10];
                n[1] = nnodes[11];
                splittingEdges.push_back(m->CreateEdge(n).first);

                auto res = Face::SplitFace(iface->getAsFace(), splittingEdges, 0);
                printf("Splitted into %zu\n", res.size());
            }
        }
    }

    m->Save("refined.vtk");
}

void createNodesForRefinedCell(Mesh *m, Mesh *mesh, Cell cell, ElementArray<Node> & nds)
{
    ElementArray<Node> nodes = cell.getNodes();
    double cnt[3];
    //printf("Cell %d has %zu nodes\n", cell.LocalID(), nodes.size());
    if(nodes.size() != 8){
        std::cout << "Cell doesn't have 8 nodes: unsupported" << std::endl;
    }

    // Determine cell boundaries
    double maxcnt[3], mincnt[3];

    std::sort(nodes.begin(), nodes.end(), Mesh::CentroidComparator(m));

    nodes[0].Barycenter(mincnt);
    nodes[7].Barycenter(maxcnt);

    double xmin = mincnt[0];
    double ymin = mincnt[1];
    double zmin = mincnt[2];
    double xmax = maxcnt[0];
    double ymax = maxcnt[1];
    double zmax = maxcnt[2];
    double dx = (xmax - xmin)/2.;
    double dy = (ymax - ymin)/2.;
    double dz = (zmax - zmin)/2.;

    //printf("dx = %e, dy = %e, dz = %e\n", dx, dy, dz);

    //printf("Cell %d: [%lf %lf]  x  [%lf %lf]  x  [%lf %lf]\n", cell.LocalID(), xmin, xmax, ymin, ymax, zmin, zmax);

    // Determine nodes 0, 2, 6, 8, 18, 20, 24, 26
    double eps = 1e-7;
    for(auto inode = nodes.begin(); inode != nodes.end(); inode++){
        inode->Barycenter(cnt);
        if(fabs(cnt[0] - xmin) < eps && fabs(cnt[1] - ymin) < eps && fabs(cnt[2] - zmin) < eps){
            inode->Barycenter(cnt);
            nds[0] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmax) < eps && fabs(cnt[1] - ymin) < eps && fabs(cnt[2] - zmin) < eps){
            inode->Barycenter(cnt);
            nds[2] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmin) < eps && fabs(cnt[1] - ymax) < eps && fabs(cnt[2] - zmin) < eps){
            inode->Barycenter(cnt);
            nds[6] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmax) < eps && fabs(cnt[1] - ymax) < eps && fabs(cnt[2] - zmin) < eps){
            inode->Barycenter(cnt);
            nds[8] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmin) < eps && fabs(cnt[1] - ymin) < eps && fabs(cnt[2] - zmax) < eps){
            inode->Barycenter(cnt);
            nds[18] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmax) < eps && fabs(cnt[1] - ymin) < eps && fabs(cnt[2] - zmax) < eps){
            inode->Barycenter(cnt);
            nds[20] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmin) < eps && fabs(cnt[1] - ymax) < eps && fabs(cnt[2] - zmax) < eps){
            inode->Barycenter(cnt);
            nds[24] = mesh->CreateNode(cnt);
        }
        if(fabs(cnt[0] - xmax) < eps && fabs(cnt[1] - ymax) < eps && fabs(cnt[2] - zmax) < eps){
            inode->Barycenter(cnt);
            nds[26] = mesh->CreateNode(cnt);
        }
    }

    // Manual specification for the rest
    {

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin;
        nds[1] = mesh->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nds[3] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nds[4] = mesh->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin;
        nds[5] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin;
        nds[7] = mesh->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nds[9] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nds[10] = mesh->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin;
        cnt[2] = zmin + dz;
        nds[11] = mesh->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nds[12] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nds[13] = mesh->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmin + dz;
        nds[14] = mesh->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nds[15] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nds[16] = mesh->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymax;
        cnt[2] = zmin + dz;
        nds[17] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin;
        cnt[2] = zmax;
        nds[19] = mesh->CreateNode(cnt);

        cnt[0] = xmin;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nds[21] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nds[22] = mesh->CreateNode(cnt);

        cnt[0] = xmax;
        cnt[1] = ymin + dy;
        cnt[2] = zmax;
        nds[23] = mesh->CreateNode(cnt);

        cnt[0] = xmin + dx;
        cnt[1] = ymax;
        cnt[2] = zmax;
        nds[25] = mesh->CreateNode(cnt);

        //std::cout << "Created nodes for cell " << cell.LocalID() << std::endl;
    }
}

void Refine2(Mesh *m)
{
    Mesh *mesh = new Mesh();
    double cnt[3];

    std::vector<Tag> tags, tagsold;
    for(unsigned i = 0; i < tagNames.size(); i++){
        if(!m->HaveTag(tagNames[i])){
            std::cout << "Input mesh has no tag '" << tagNames[i] << "', skipping" << std::endl;
            continue;
        }
        tagsold.push_back(m->GetTag(tagNames[i]));
        tags.push_back(mesh->CreateTag(tagNames[i], DATA_REAL, CELL, NONE, 1));
        std::cout << "Added tag '" << tagNames[i] << "'\n";
    }
    unsigned nT = tags.size();

    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        //printf("Processing cell %d\n", icell->LocalID());
        ElementArray<Node> nnodes(mesh, 27);
        createNodesForRefinedCell(m, mesh, icell->getAsCell(), nnodes);

        ElementArray<Face> newFaces(mesh);
        ElementArray<Node> faceNodes(mesh, 4);

        newFaces.clear();
        Face f;
        faceNodes[0] = nnodes[0];
        faceNodes[1] = nnodes[1];
        faceNodes[2] = nnodes[10];
        faceNodes[3] = nnodes[9];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        faceNodes[0] = nnodes[1];
        faceNodes[1] = nnodes[4];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[10];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        faceNodes[0] = nnodes[4];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[12];
        faceNodes[3] = nnodes[3];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        faceNodes[0] = nnodes[3];
        faceNodes[1] = nnodes[12];
        faceNodes[2] = nnodes[9];
        faceNodes[3] = nnodes[0];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        faceNodes[0] = nnodes[0];
        faceNodes[1] = nnodes[1];
        faceNodes[2] = nnodes[4];
        faceNodes[3] = nnodes[3];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        faceNodes[0] = nnodes[9];
        faceNodes[1] = nnodes[10];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[12];
        f = mesh->CreateFace(faceNodes).first;
        //f.Barycenter(cnt);
        //printf("f: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        newFaces.push_back(f);
        Cell c = mesh->CreateCell(newFaces).first;
        //c.Barycenter(cnt);
        //printf("Cell cnt: %e %e %e\n", cnt[0], cnt[1], cnt[2]);
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[1];
        faceNodes[1] = nnodes[2];
        faceNodes[2] = nnodes[11];
        faceNodes[3] = nnodes[10];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[2];
        faceNodes[1] = nnodes[5];
        faceNodes[2] = nnodes[14];
        faceNodes[3] = nnodes[11];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[5];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[4];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[4];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[10];
        faceNodes[3] = nnodes[1];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[1];
        faceNodes[1] = nnodes[2];
        faceNodes[2] = nnodes[5];
        faceNodes[3] = nnodes[4];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[10];
        faceNodes[1] = nnodes[11];
        faceNodes[2] = nnodes[14];
        faceNodes[3] = nnodes[13];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[3];
        faceNodes[1] = nnodes[4];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[12];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[4];
        faceNodes[1] = nnodes[7];
        faceNodes[2] = nnodes[16];
        faceNodes[3] = nnodes[13];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[7];
        faceNodes[1] = nnodes[16];
        faceNodes[2] = nnodes[15];
        faceNodes[3] = nnodes[6];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[6];
        faceNodes[1] = nnodes[15];
        faceNodes[2] = nnodes[12];
        faceNodes[3] = nnodes[3];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[3];
        faceNodes[1] = nnodes[4];
        faceNodes[2] = nnodes[7];
        faceNodes[3] = nnodes[6];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[12];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[16];
        faceNodes[3] = nnodes[15];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[4];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[14];
        faceNodes[3] = nnodes[5];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[5];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[17];
        faceNodes[3] = nnodes[8];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[8];
        faceNodes[1] = nnodes[17];
        faceNodes[2] = nnodes[16];
        faceNodes[3] = nnodes[7];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[7];
        faceNodes[1] = nnodes[16];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[4];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[4];
        faceNodes[1] = nnodes[5];
        faceNodes[2] = nnodes[8];
        faceNodes[3] = nnodes[7];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[17];
        faceNodes[3] = nnodes[16];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[9];
        faceNodes[1] = nnodes[10];
        faceNodes[2] = nnodes[19];
        faceNodes[3] = nnodes[18];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[10];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[22];
        faceNodes[3] = nnodes[19];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[22];
        faceNodes[2] = nnodes[21];
        faceNodes[3] = nnodes[12];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[12];
        faceNodes[1] = nnodes[21];
        faceNodes[2] = nnodes[18];
        faceNodes[3] = nnodes[9];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[9];
        faceNodes[1] = nnodes[10];
        faceNodes[2] = nnodes[13];
        faceNodes[3] = nnodes[12];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[18];
        faceNodes[1] = nnodes[19];
        faceNodes[2] = nnodes[22];
        faceNodes[3] = nnodes[21];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[10];
        faceNodes[1] = nnodes[11];
        faceNodes[2] = nnodes[20];
        faceNodes[3] = nnodes[19];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[11];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[23];
        faceNodes[3] = nnodes[20];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[14];
        faceNodes[1] = nnodes[23];
        faceNodes[2] = nnodes[22];
        faceNodes[3] = nnodes[13];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[22];
        faceNodes[2] = nnodes[19];
        faceNodes[3] = nnodes[10];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[10];
        faceNodes[1] = nnodes[11];
        faceNodes[2] = nnodes[14];
        faceNodes[3] = nnodes[13];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[19];
        faceNodes[1] = nnodes[20];
        faceNodes[2] = nnodes[23];
        faceNodes[3] = nnodes[22];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[12];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[22];
        faceNodes[3] = nnodes[21];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[16];
        faceNodes[2] = nnodes[25];
        faceNodes[3] = nnodes[22];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[16];
        faceNodes[1] = nnodes[25];
        faceNodes[2] = nnodes[24];
        faceNodes[3] = nnodes[15];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[15];
        faceNodes[1] = nnodes[24];
        faceNodes[2] = nnodes[21];
        faceNodes[3] = nnodes[12];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[12];
        faceNodes[1] = nnodes[13];
        faceNodes[2] = nnodes[16];
        faceNodes[3] = nnodes[15];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[21];
        faceNodes[1] = nnodes[22];
        faceNodes[2] = nnodes[25];
        faceNodes[3] = nnodes[24];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }


        newFaces.clear();
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[23];
        faceNodes[3] = nnodes[22];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[14];
        faceNodes[1] = nnodes[17];
        faceNodes[2] = nnodes[26];
        faceNodes[3] = nnodes[23];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[17];
        faceNodes[1] = nnodes[26];
        faceNodes[2] = nnodes[25];
        faceNodes[3] = nnodes[16];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[16];
        faceNodes[1] = nnodes[25];
        faceNodes[2] = nnodes[22];
        faceNodes[3] = nnodes[13];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[13];
        faceNodes[1] = nnodes[14];
        faceNodes[2] = nnodes[17];
        faceNodes[3] = nnodes[16];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        faceNodes[0] = nnodes[22];
        faceNodes[1] = nnodes[23];
        faceNodes[2] = nnodes[26];
        faceNodes[3] = nnodes[25];
        f = mesh->CreateFace(faceNodes).first;
        newFaces.push_back(f);
        c = mesh->CreateCell(newFaces).first;
        for(unsigned k = 0; k < nT; k++){
            c.Real(tags[k]) = icell->Real(tagsold[k]);
        }
    }
    printf("Mesh has %d cells\n", mesh->NumberOfCells());
    mesh->AssignGlobalID(CELL);
    mesh->Save("refined.vtk");
    delete mesh;
}

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: gridrefine <mesh_file> <refinement_level>" << std::endl;
        exit(1);
    }
    std::cout << "Note: Only cubic-like grids are supported" << std::endl;

    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);

    //Mesh *mesh = new Mesh("mesh");
    Mesh *m    = new Mesh("refined mesh");

    /*int rank = mesh->GetProcessorRank();
    std::cout << "rank = " << rank << std::endl;
    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    if(rank == 0)
        mesh->Load(argv[1]);

    Partitioner p(mesh);
    p.SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
    p.Evaluate();

    mesh->Redistribute();
    mesh->ExchangeGhost(1, NODE);

    Tag T = mesh->CreateTag("My", DATA_REAL, CELL, false, 1);
    T.SetPrint(true);
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost)
            continue;
        double x[3];
        icell->Barycenter(x);
        icell->Real(T) = x[0]+x[1]+x[2];
    }*/
    m->Load(argv[1]);

    tagNames.push_back("K");
    tagNames.push_back("Permeability_scalar");
    tagNames.push_back("PORO");

    Refine2(m);

    delete m;

    Mesh::Finalize();
    Partitioner::Finalize();
}
