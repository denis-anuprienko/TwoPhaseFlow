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

void ReduceMax(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
{
    (void) size;
    element->Integer(tag) = std::max(element->Integer(tag),*((const INMOST_DATA_INTEGER_TYPE *)data));
}

void RefineA(Mesh *m)
{
    std::vector<std::string> tagNames;
    tagNames.push_back("PORO");
    tagNames.push_back("Permeability_scalar");
    tagNames.push_back("Permeability");

    TagInteger indicator = m->CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);

    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++)
        icell->Integer(indicator) = 1;

    std::vector<Tag> tags;
    for(unsigned i = 0; i < tagNames.size(); i++){
        if(!m->HaveTag(tagNames[i])){
            std::cout << "Input mesh has no tag '" << tagNames[i] << "', skipping" << std::endl;
            continue;
        }
        tags.push_back(m->GetTag(tagNames[i]));
        //tags.push_back(mesh->CreateTag(tagNames[i], DATA_REAL, CELL, NONE, 1));
        std::cout << "Added tag '" << tagNames[i] << "'\n";
    }
    int nT = tags.size();

    static int call_counter = 0;
    int ret = 0; //return number of refined cells

    ////model = NULL;
    //create a tag that stores maximal refinement level of each element
    TagInteger level = m->CreateTag("REFINEMENT_LEVEL",DATA_INTEGER,CELL|FACE|EDGE|NODE|ESET,NONE,1);
    //tag_status = m->CreateTag("TAG_STATUS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
    TagInteger set_id = m->CreateTag("SET_ID",DATA_INTEGER,CELL|ESET,ESET,1);
    //tag_an = m->CreateTag("TAG_AN",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
    //ref_tag = m->CreateTag("REF",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
    //create a tag that stores links to all the hanging nodes of the cell
    TagReferenceArray hanging_nodes = m->CreateTag("HANGING_NODES",DATA_REFERENCE,CELL|FACE,NONE);
    //create a tag that stores links to sets
    TagReference parent_set = m->CreateTag("PARENT_SET",DATA_REFERENCE,CELL,NONE,1);
    int size = m->GetProcessorsNumber();
    int rank = m->GetProcessorRank();

    ElementSet root;
    if( !root.isValid() )
    {
        root = m->GetSet("AM_ROOT_SET");
        if( root == InvalidElement() )
        {
            root = m->CreateSetUnique("AM_ROOT_SET").first;
            //root.SetExchange(ElementSet::SYNC_ELEMENTS_SHARED);
            level[root] = 0;
            for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
            {
                root.PutElement(it->self());
                parent_set[it->self()] = root.GetHandle();
            }
            m->Enumerate(CELL,set_id);
        }
    }
    if( !m->HaveGlobalID(CELL) ) m->AssignGlobalID(CELL);

    //m->CheckSetLinks(__FILE__,__LINE__);
    //CheckParentSet(__FILE__,__LINE__);


    int schedule_counter = 1; //indicates order in which refinement will be scheduled
    int scheduled = 1; //indicates that at least one element was scheduled on current sweep
    //0. Extend indicator for edges and faces
    indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
    //1.Communicate indicator - it may be not synced
    m->ExchangeData(indicator,CELL,0);
    while(scheduled)
    {
        //2.Propogate indicator down to the faces,edges
        //  select schedule for them
        for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
        {
            Cell c = m->CellByLocalID(it);
            if( indicator[c] == schedule_counter )
            {
                ElementArray<Element> adj = c.getAdjElements(FACE|EDGE);
                for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
                {
                    if( level[adj[kt]] == level[c] ) //do not schedule finer or coarser elements
                        indicator[adj[kt]] = schedule_counter; //refine together with current cell
                }
            }
        }
        //3.Communicate indicator on faces and edges
        ////m->ReduceData(indicator,FACE|EDGE,0,ReduceMax);
        m->ExchangeData(indicator,FACE|EDGE,0);
        //4.Check for each cell if there is
        //  any hanging node with adjacent in a need to refine,
        //  schedule for refinement earlier.
        scheduled = 0;
        for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
        {
            Cell c = m->CellByLocalID(it);
            //already scheduled cells may be required to be refined first
            //if( indicator[c] == 0 ) //some optimization
            {
                bool scheduled_c = false;
                //any finer level edge is scheduled to be refined first
                ElementArray<Edge> edges = c->getEdges();
                for(ElementArray<Edge>::size_type kt = 0; kt < edges.size() && !scheduled_c; ++kt)
                {
                    //if a finer edge is scheduled
                    //then this cell should be refined first
                    if( indicator[edges[kt]] != 0 &&
                        level[edges[kt]] > level[c] &&
                        indicator[edges[kt]] >= indicator[c] )
                    {
                        indicator[c] = schedule_counter+1;
                        scheduled++;
                        scheduled_c = true;
                    }
                }
            }
        }
        //5.Exchange indicator on cells
        m->ReduceData(indicator,CELL,0,ReduceMax);
        m->ExchangeData(indicator,CELL,0);
        //6.Go back to 1 until no new elements scheduled
        scheduled = m->Integrate(scheduled);
        if( scheduled ) schedule_counter++;
    }
    //m->ExchangeData(indicator,CELL | FACE | EDGE,0);
    printf("First while loop finished\n");


    m->ExchangeData(hanging_nodes,CELL | FACE,0);


    //6.Refine
    m->BeginModification();
    printf("Began mod\n");
    int iti = 0;
    while(schedule_counter)
    {
        printf("schedule_counter loop it %d\n", iti);
        iti++;

        Storage::real xyz[3] = {0,0,0};
        //7.split all edges of the current schedule
        for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
        {
            Edge e = m->EdgeByLocalID(it);
            if( !e.Hidden() && indicator[e] == schedule_counter )
            {
                //remember adjacent faces that should get information about new hanging node
                ElementArray<Face> edge_faces = e.getFaces();
                //location on the center of the edge
                for(Storage::integer d = 0; d < m->GetDimensions(); ++d)
                    xyz[d] = (e.getBeg().Coords()[d]+e.getEnd().Coords()[d])*0.5;

                // create middle node
                Node n = m->CreateNode(xyz);

                //set increased level for new node
                level[n] = level[e.getBeg()] = level[e.getEnd()] = level[e]+1;

                //for each face provide link to a new hanging node
                for(ElementArray<Face>::size_type kt = 0; kt < edge_faces.size(); ++kt)
                    hanging_nodes[edge_faces[kt]].push_back(n);

                //split the edge by the middle node
                ElementArray<Edge> new_edges = Edge::SplitEdge(e,ElementArray<Node>(m,1,n.GetHandle()),0);

                //set increased level for new edges
                level[new_edges[0]] = level[new_edges[1]] = level[e]+1;
            }
        }
        printf("  Processed all edges\n");

        //8.split all faces of the current schedule, using hanging nodes on edges
        for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
        {
            Face f = m->FaceByLocalID(it);
            if( !f.Hidden() && indicator[f] == schedule_counter )
            {
                //connect face center to hanging nodes of the face
                Storage::reference_array face_hanging_nodes = hanging_nodes[f];
                //remember adjacent cells that should get information about new hanging node
                //and new hanging edges
                ElementArray<Cell> face_cells = f.getCells();
                //create node at face center
                //f->Centroid(xyz);
                for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
                for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
                    for(int d = 0; d < 3; ++d) xyz[d] += face_hanging_nodes[kt].getAsNode().Coords()[d];
                for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)face_hanging_nodes.size();
                //todo: request transformation of node location according to geometrical model
                //create middle node
                Node n = m->CreateNode(xyz);
                //set increased level for the new node
                level[n] = level[f]+1;
                //for each cell provide link to new hanging node
                for(ElementArray<Face>::size_type kt = 0; kt < face_cells.size(); ++kt)
                    hanging_nodes[face_cells[kt]].push_back(n);
                ElementArray<Node> edge_nodes(m,2); //to create new edges
                ElementArray<Edge> hanging_edges(m,face_hanging_nodes.size());
                edge_nodes[0] = n;
                for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
                {
                    edge_nodes[1] = face_hanging_nodes[kt].getAsNode();
                    hanging_edges[kt] = m->CreateEdge(edge_nodes).first;
                    //set increased level for new edges
                    level[hanging_edges[kt]] = level[f]+1;
                }
                //split the face by these edges
                ElementArray<Face> new_faces = Face::SplitFace(f,hanging_edges,0);
                //set increased level to new faces
                for(ElementArray<Face>::size_type kt = 0; kt < new_faces.size(); ++kt)
                    level[new_faces[kt]] = level[f]+1;
            }
        }

        printf("  Processed faces\n");
        //this tag helps recreate internal face
        TagReferenceArray internal_face_edges = m->CreateTag("INTERNAL_FACE_EDGES",DATA_REFERENCE,NODE,NODE,4);
        //this marker helps detect edges of current cell only
        MarkerType mark_cell_edges = m->CreateMarker();
        //this marker helps detect nodes hanging on edges of unrefined cell
        MarkerType mark_hanging_nodes = m->CreateMarker();
        //9.split all cells of the current schedule
        for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
        {
            Cell c = m->CellByLocalID(it);

            //printf("Processing cell %d\n", it);
            if( !c.Hidden() && indicator[c] == schedule_counter )
            {
                Storage::reference_array cell_hanging_nodes = hanging_nodes[c]; //nodes to be connected
                //create node at cell center
                for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
                for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
                    for(int d = 0; d < 3; ++d) xyz[d] += cell_hanging_nodes[kt].getAsNode().Coords()[d];
                for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)cell_hanging_nodes.size();


                //create middle node
                Node n = m->CreateNode(xyz);

                //set increased level for the new node
                level[n] = level[c]+1;

                //retrive all edges of current face to mark them
                ElementArray<Edge> cell_edges = c.getEdges();
#if !defined(NDEBUG)
                for(ElementArray<Edge>::iterator jt = cell_edges.begin(); jt != cell_edges.end(); ++jt) assert(level[*jt] == level[c]+1);
                ElementArray<Face> cell_faces = c.getFaces();
                for(ElementArray<Face>::iterator jt = cell_faces.begin(); jt != cell_faces.end(); ++jt) assert(level[*jt] == level[c]+1);
#endif //NDEBUG
                //mark all edges so that we can retive them later
                cell_edges.SetMarker(mark_cell_edges);

                //connect face center to centers of faces by edges
                ElementArray<Node> edge_nodes(m,2);
                ElementArray<Edge> edges_to_faces(m,cell_hanging_nodes.size());
                edge_nodes[0] = n;
                for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
                {
                    assert(cell_hanging_nodes[kt].isValid());
                    //todo: unmark hanging node on edge if no more cells depend on it
                    edge_nodes[1] = cell_hanging_nodes[kt].getAsNode();
                    edges_to_faces[kt] = m->CreateEdge(edge_nodes).first;
                    //set increased level for new edges
                    level[edges_to_faces[kt]] = level[c]+1;
                    //for each node other then the hanging node of the face
                    //(this is hanging node on the edge)
                    //we record a pair of edges to reconstruct internal faces
                    ElementArray<Edge> hanging_edges = cell_hanging_nodes[kt].getEdges(mark_cell_edges,0);
                    for(ElementArray<Edge>::size_type lt = 0; lt < hanging_edges.size(); ++lt)
                    {
                        //get hanging node on the edge
                        assert(hanging_edges[lt].getBeg() == cell_hanging_nodes[kt] || hanging_edges[lt].getEnd() == cell_hanging_nodes[kt]);
                        Node v = hanging_edges[lt].getBeg() == cell_hanging_nodes[kt]? hanging_edges[lt].getEnd() : hanging_edges[lt].getBeg();
                        //mark so that we can collect all of them
                        v.SetMarker(mark_hanging_nodes);
                        //fill the edges
                        Storage::reference_array face_edges = internal_face_edges[v];
                        //fill first two in forward order
                        //this way we make a closed loop
                        assert(face_edges[0] == InvalidElement() || face_edges[2] == InvalidElement());
                        if( face_edges[0] == InvalidElement() )
                        {
                            face_edges[0] = edges_to_faces[kt];
                            face_edges[1] = hanging_edges[lt];
                        }
                        else //last two in reverse
                        {
                            assert(face_edges[2] ==InvalidElement());
                            face_edges[2] = hanging_edges[lt];
                            face_edges[3] = edges_to_faces[kt];
                        }
                    }
                }

                //printf("Connected face centers\n");
                //remove marker from cell edges
                cell_edges.RemMarker(mark_cell_edges);
                //now we have to create internal faces
                ElementArray<Node> edge_hanging_nodes = c.getNodes(mark_hanging_nodes,0);
                ElementArray<Face> internal_faces(m,edge_hanging_nodes.size());
                //unmark hanging nodes on edges
                edge_hanging_nodes.RemMarker(mark_hanging_nodes);
                for(ElementArray<Node>::size_type kt = 0; kt < edge_hanging_nodes.size(); ++kt)
                {
                    //create a face based on collected edges
                    Storage::reference_array face_edges = internal_face_edges[edge_hanging_nodes[kt]];
                    assert(face_edges[0].isValid());
                    assert(face_edges[1].isValid());
                    assert(face_edges[2].isValid());
                    assert(face_edges[3].isValid());
                    internal_faces[kt] = m->CreateFace(ElementArray<Edge>(m,face_edges.begin(),face_edges.end())).first;
                    //set increased level
                    level[internal_faces[kt]] = level[c]+1;
                    //clean up structure, so that other cells can use it
                    edge_hanging_nodes[kt].DelData(internal_face_edges);
                }

                //split the cell
                //retrive parent set
                ElementSet parent(m,parent_set[c]);
                //create set corresponding to old coarse cell
                //Storage::real cnt[3];
                //c.Centroid(cnt);
                std::stringstream set_name;
                //set_name << parent.GetName() << "_C" << c.GlobalID(); //rand may be unsafe
                if( parent == root )
                    set_name << "AM_R" << set_id[c];
                else
                    set_name << parent.GetName() << "C" << set_id[c];
                //set_name << base64_encode_((unsigned char *)cnt,3*sizeof(double)/sizeof(unsigned char));

                ElementSet check_set = m->GetSet(set_name.str());
                if( check_set.isValid() )
                {
                    std::cout << rank << " set " << set_name.str() << " for cell " << c.GlobalID() << " " << Element::StatusName(c.GetStatus()) << " already exists" << std::endl;
                    if( check_set->HaveParent() )
                        std::cout << rank << " parent is " << check_set->GetParent()->GetName() << " cell parent is " << parent.GetName() << std::endl;
                    std::cout << rank << " Elements of " << check_set.GetName() << ": ";
                    for(ElementSet::iterator it = check_set.Begin(); it != check_set.End(); ++it)
                        std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << "," << indicator[*it] << " ";
                    std::cout << std::endl;
                    std::cout << rank << " Elements of " << parent.GetName() << ": ";
                    for(ElementSet::iterator it = parent.Begin(); it != parent.End(); ++it)
                        std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << "," << indicator[*it] << " ";
                    std::cout << std::endl;
                    if( parent.HaveChild() )
                    {
                        std::cout << rank << " Children of " << parent.GetName() << ": ";
                        for(ElementSet jt = parent.GetChild(); jt.isValid(); jt = jt.GetSibling() )
                            std::cout << jt.GetName() << " size " << jt.Size() << " ";
                        std::cout << std::endl;
                    }
                    exit(-1);
                }

                ElementSet cell_set = m->CreateSetUnique(set_name.str()).first;
                //cell_set->SetExchange(ElementSet::SYNC_ELEMENTS_ALL);
                level[cell_set] = level[c]+1;
                set_id[cell_set] = set_id[c];

                ElementArray<Cell> new_cells = Cell::SplitCell(c,internal_faces,0);
                std::sort(new_cells.begin(),new_cells.end(),Mesh::CentroidComparator(m));

                for(auto icc = new_cells.begin(); icc != new_cells.end(); icc++){
                    for(int k = 0; k < nT; k++){
                        icc->Real(tags[k]) = c.Real(tags[k]);
                    }
                }

                //set up increased level for the new cells
                for(ElementArray<Cell>::size_type kt = 0; kt < new_cells.size(); ++kt)
                {
                    set_id[new_cells[kt]] = kt;
                    level[new_cells[kt]] = level[c]+1;
                    cell_set.PutElement(new_cells[kt]);
                    parent_set[new_cells[kt]] = cell_set.GetHandle();
                }

                //parent.AddChild(cell_set);
                //printf("Added ch\n");
                ret++;
            }
        }

        printf("Processed all cells\n");
        m->ReleaseMarker(mark_hanging_nodes);
        m->ReleaseMarker(mark_cell_edges);
        m->DeleteTag(internal_face_edges);
        //10.jump to later schedule, and go to 7.
        schedule_counter--;
    }
    m->CheckSetLinks(__FILE__,__LINE__);
    //free created tag
    m->DeleteTag(indicator);
    m->DeleteTag(set_id);
    m->DeleteTag(parent_set);
    m->DeleteTag(level);


    m->CheckSetLinks(__FILE__,__LINE__);

    //11. Restore parallel connectivity, global ids
    m->ResolveModification();

    //12. Let the user update their data
    //if( model ) model->Adaptation(*m);
    m->CheckSetLinks(__FILE__,__LINE__);
    //13. Delete old elements of the mesh
    m->ApplyModification();

    //m->ExchangeGhost(1,NODE,m->NewMarker());
    //14. Done
    //cout << rank << ": Before end " << std::endl;
    m->EndModification();

    //keep links to prevent loss during balancing
    m->ExchangeData(parent_set,CELL,0);
    m->ExchangeData(hanging_nodes,CELL | FACE,0);

    ////CheckParentSet(__FILE__,__LINE__);

    //reorder element's data to free up space
    m->ReorderEmpty(CELL|FACE|EDGE|NODE);
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

    m->Load(argv[1]);

    tagNames.push_back("K");
    tagNames.push_back("Permeability_scalar");
    tagNames.push_back("Permeability");
    tagNames.push_back("PORO");

    //Refine2(m);

    RefineA(m);

    std::cout << "Generated mesh with " << m->NumberOfCells() << " cells" << std::endl;

    m->Save("seed.vtk");

    delete m;

    Mesh::Finalize();
    Partitioner::Finalize();
}
