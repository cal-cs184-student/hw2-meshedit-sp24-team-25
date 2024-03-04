#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
  Vector2D lerp(Vector2D a, Vector2D b, float t) {
    return a + (b - a) * t;
  }

  Vector3D lerp(Vector3D a, Vector3D b, float t) {
    return a + (b - a) * t;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    std::vector<Vector2D> output;
    for (int i = 0; i < points.size() - 1; i++) {
      output.push_back(lerp(points[i], points[i + 1], t));
    }
    return output;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> output;
    for (int i = 0; i < points.size() - 1; i++) {
      output.push_back(lerp(points[i], points[i + 1], t));
    }
    return output;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    if (points.size() == 1) {
      return points[0];
    }
    return evaluate1D(evaluateStep(points, t), t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    std::vector<Vector3D> col;
    for (std::vector<Vector3D> row : controlPoints) {
      col.push_back(evaluate1D(row, u));
    }
    return evaluate1D(col, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    Vector3D normal;
    
    HalfedgeCIter h = halfedge();      // get the outgoing half-edge of the vertex
    do {
      HalfedgeCIter e2 = h->next();
      VertexCIter v2 = e2->vertex();
      
      HalfedgeCIter e3 = e2->next();
      VertexCIter v3 = e3->vertex();
      
      Vector3D ab = v2->position - position;
      Vector3D ac = v3->position - position;
      
      
      normal += cross(ab, ac);
      
      HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
      h = h_twin->next();               // move to the next outgoing half-edge of the vertex
    } while(h != halfedge());          // keep going until we are back where we were
    
    normal.normalize();
    return normal;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    
    if (e0->isBoundary()) {
      return e0;
    }
    
    HalfedgeIter bc = e0->halfedge();
    HalfedgeIter ca = bc->next();
    HalfedgeIter ab = ca->next();
    
    HalfedgeIter cb = bc->twin();
    HalfedgeIter bd = cb->next();
    HalfedgeIter dc = bd->next();
    
    VertexIter a = ab->vertex();
    VertexIter b = bc->vertex();
    VertexIter c = ca->vertex();
    VertexIter d = dc->vertex();
    
    FaceIter abc = ab->face();
    FaceIter bdc = bd->face();
    
    a->halfedge() = ab;
    b->halfedge() = bd;
    c->halfedge() = ca;
    d->halfedge() = dc;
    
    abc->halfedge() = ca;
    bdc->halfedge() = bd;
    
    ab->setNeighbors(bd, ab->twin(), a, ab->edge(), bdc);
    bd->setNeighbors(cb, bd->twin(), b, bd->edge(), bdc);
    cb->setNeighbors(ab, bc, d, e0, bdc);
    
    dc->setNeighbors(ca, dc->twin(), d, dc->edge(), abc);
    ca->setNeighbors(bc, ca->twin(), c, ca->edge(), abc);
    bc->setNeighbors(dc, cb, a, e0, abc);
    
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    
    if (e0->isBoundary()) {
      return VertexIter();
    }
    
    HalfedgeIter bc = e0->halfedge();
    HalfedgeIter ca = bc->next();
    HalfedgeIter ab = ca->next();
    
    HalfedgeIter cb = bc->twin();
    HalfedgeIter bd = cb->next();
    HalfedgeIter dc = bd->next();
    
    VertexIter a = ab->vertex();
    VertexIter b = bc->vertex();
    VertexIter c = ca->vertex();
    VertexIter d = dc->vertex();
    
//    EdgeIter e_bc = bc->edge();
//    EdgeIter e_ca = ca->edge();
//    EdgeIter e_ab = ab->edge();
//    
//    EdgeIter e_cb = cb->edge();
//    EdgeIter e_bd = bd->edge();
//    EdgeIter e_dc = dc->edge();
    
    FaceIter abc = ab->face();
    FaceIter bdc = bd->face();
    
    VertexIter m = newVertex();
  
    HalfedgeIter am = newHalfedge();
    HalfedgeIter ma = newHalfedge();
    HalfedgeIter bm = newHalfedge();
    HalfedgeIter mb = newHalfedge();
    HalfedgeIter dm = newHalfedge();
    HalfedgeIter md = newHalfedge();
    
    EdgeIter e_am = newEdge();
    EdgeIter e_bm = newEdge();
    EdgeIter e_dm = newEdge();
    
    FaceIter abm = newFace();
    FaceIter bdm = newFace();
    
    a->halfedge() = ab;
    b->halfedge() = bd;
    c->halfedge() = ca;
    d->halfedge() = dc;
    
    m->position = (c->position + b->position) / 2;
    m->halfedge() = bc;
    m->isNew = true;
    m->newPosition = e0->newPosition;
    
    bc->setNeighbors(ca, cb, m, bc->edge(), abc);
    ca->setNeighbors(am, ca->twin(), c, ca->edge(), abc);
    am->setNeighbors(bc, ma, a, e_am, abc);
    
    ab->setNeighbors(bm, ab->twin(), a, ab->edge(), abm);
    bm->setNeighbors(ma, mb, b, e_bm, abm);
    ma->setNeighbors(ab, am, m, e_am, abm);
    
    bd->setNeighbors(dm, bd->twin(), b, bd->edge(), bdm);
    dm->setNeighbors(mb, md, d, e_dm, bdm);
    mb->setNeighbors(bd, bm, m, e_bm, bdm);
    
    dc->setNeighbors(cb, dc->twin(), d, dc->edge(), bdc);
    cb->setNeighbors(md, bc, c, cb->edge(), bdc);
    md->setNeighbors(dc, dm, m, e_dm, bdc);
    
    e_am->halfedge() = am;
    e_bm->halfedge() = bm;
    e_dm->halfedge() = dm;
    
    e_am->isSplit = true;
    e_bm->isSplit = true;
    e_dm->isSplit = true;
    bc->edge()->isSplit = true;
    
    e_am->isNew = true;
    e_dm->isNew = true;
    
    abc->halfedge() = ca;
    abm->halfedge() = ab;
    bdm->halfedge() = bd;
    bdc->halfedge() = dc;
    
    return VertexIter();
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
{
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.
    
    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.
    
    // 5. Copy the new vertex positions into final Vertex::position.
    
    VertexIter v = mesh.verticesBegin();
    while (v != mesh.verticesEnd()) {
        
      // get the next edge NOW!
      VertexIter nextVertex = v;
      nextVertex++;
      
      HalfedgeIter h = v->halfedge();
        
      Vector3D position;
      int n = 0;
      do {
        HalfedgeCIter e2 = h->next();
        VertexCIter v2 = e2->vertex();
        position += v2->position;
        
        n += 1;
        HalfedgeIter h_twin = h->twin();   // get the opposite half-edge
        h = h_twin->next();                 // move to the next outgoing half-edge of the vertex
      } while(h != v->halfedge());          // keep going until we are back where we were
      
      float u;
      if (n == 3) {
        u = 3.0/16;
      } else {
        u = 3 / (8.0 * n);
      }
      
//      cout << u << endl;
      
      position *= u;
      
      position += (1 - n * u) * v->position;
      v->newPosition = position;
      
      // ... we still have a valid reference to the next edge
      v = nextVertex;
    }
    
    
    EdgeIter e = mesh.edgesBegin();
    while (e != mesh.edgesEnd()) {
      
      // get the next edge NOW!
      EdgeIter nextEdge = e;
      nextEdge++;
      
      Vector3D edgePos;
      HalfedgeCIter h = e->halfedge();
      
      edgePos += h->vertex()->position * 3.0/8;
      h = h->next();
      edgePos += h->vertex()->position * 3.0/8;
      h = h->next();
      edgePos += h->vertex()->position * 1.0/8;
      
      h = e->halfedge();
      h = h->twin()->next()->next();
      edgePos += h->vertex()->position * 1.0/8;
      
      e->newPosition = edgePos;
      
      // ... we still have a valid reference to the next edge
      e = nextEdge;
    }
    
    e = mesh.edgesBegin();
    while (e != mesh.edgesEnd()) {
      // get the next edge NOW!
      EdgeIter nextEdge = e;
      nextEdge++;
      
      // now, even if splitting the edge deletes it ...
      if (!e->isSplit) {
          mesh.splitEdge(e);
      }
      
      // ... we still have a valid reference to the next edge
      e = nextEdge;
    }
    
    e = mesh.edgesBegin();
    while (e != mesh.edgesEnd()) {
      // get the next edge NOW!
      EdgeIter nextEdge = e;
      nextEdge++;
      
      // now, even if splitting the edge deletes it ...
      if (e->isNew) {
        HalfedgeIter h = e->halfedge();
        if (h->vertex()->isNew != h->twin()->vertex()->isNew) {
          mesh.flipEdge(e);
        }
      }
      
      e->isNew = false;
      e->isSplit = false;
      
      // ... we still have a valid reference to the next edge
      e = nextEdge;
    }
    
    v = mesh.verticesBegin();
    while (v != mesh.verticesEnd()) {
        
      // get the next edge NOW!
      VertexIter nextVertex = v;
      nextVertex++;
      
      v->isNew = false;
      v->position = v->newPosition;
      
      // ... we still have a valid reference to the next edge
      v = nextVertex;
    }
  }
    
    
}
