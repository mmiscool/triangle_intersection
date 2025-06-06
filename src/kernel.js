
/**
 * A class representing a triangle in 3D space, defined by three points {x, y, z}.
 * Provides a method to compute the intersection (if any) between two triangles.
 *
 * Intersection rules:
 * - If the triangles do not intersect, returns [].
 * - If they touch at exactly one point, returns [point].
 * - If they intersect along a line segment, returns [point1, point2].
 * - If they are coplanar and overlap, returns an array of the overlapping region’s vertices (up to 6 points),
 *   although in many typical coplanar‐triangle‐intersection cases you will see 1, 2, or 3 points returned.
 */
export class Triangle {
    /**
   * @param {{x: number, y: number, z: number}} a - first vertex
   * @param {{x: number, y: number, z: number}} b - second vertex
   * @param {{x: number, y: number, z: number}} c - third vertex
   */
    /**
   * @param {{x:number, y:number, z:number}} a - first vertex
   * @param {{x:number, y:number, z:number}} b - second vertex
   * @param {{x:number, y:number, z:number}} c - third vertex
   */
    constructor(a, b, c) {
        this.a = a;
        this.b = b;
        this.c = c;
        // Compute the (unnormalized) normal vector of the triangle’s plane: n = (b - a) × (c - a)
        this.normal = Triangle.cross(Triangle.subtract(b, a), Triangle.subtract(c, a));
        // Plane equation: n · X + d = 0  =>  d = -n · a
        this.d = -Triangle.dot(this.normal, a);
    }
    /**
   * Compute the intersection between this triangle and another triangle.
   * @param {Triangle} other - the other triangle
   * @param {number} [eps=1e-6] - epsilon tolerance for numerical comparisons
   * @returns {Array<{x: number, y: number, z: number}>} An array of 0, 1, 2, or more points.
   */
    /**
   * Compute the intersection between this triangle and another triangle.
   * @param {Triangle} other - the other triangle
   * @param {number} [eps=1e-6] - epsilon tolerance for numerical comparisons
   * @returns {Array<{x:number, y:number, z:number}>} An array of 0, 1, 2, or more points.
   */
    intersect(other, eps = 0.000001) {
        // 1) Compute signed distances of other’s vertices to this triangle’s plane
        const d1 = Triangle.dot(this.normal, other.a) + this.d;
        const d2 = Triangle.dot(this.normal, other.b) + this.d;
        const d3 = Triangle.dot(this.normal, other.c) + this.d;
        // If all three distances are > eps or all three < -eps, then other is completely on one side → no intersection
        if (d1 > eps && d2 > eps && d3 > eps)
            return [];
        if (d1 < -eps && d2 < -eps && d3 < -eps)
            return [];
        // If all distances are near zero, triangles are coplanar (within eps)
        if (Math.abs(d1) < eps && Math.abs(d2) < eps && Math.abs(d3) < eps) {
            return this._coplanarIntersection(other, eps);
        }
        // 2) Compute signed distances of this triangle’s vertices to other’s plane
        const e1 = Triangle.dot(other.normal, this.a) + other.d;
        const e2 = Triangle.dot(other.normal, this.b) + other.d;
        const e3 = Triangle.dot(other.normal, this.c) + other.d;
        if (e1 > eps && e2 > eps && e3 > eps)
            return [];
        if (e1 < -eps && e2 < -eps && e3 < -eps)
            return [];
        // 3) For non‐coplanar case, find intersection segment of each triangle with the other’s plane
        // 3a) Find points where edges of this triangle intersect other’s plane
        const triVerts1 = [
            this.a,
            this.b,
            this.c
        ];
        const distsToOtherPlane = [
            e1,
            e2,
            e3
        ];
        const ptsOnThis = this._linePlaneIntersect(other.normal, other.d, triVerts1, distsToOtherPlane, eps);
        // 3b) Find points where edges of other triangle intersect this plane
        const triVerts2 = [
            other.a,
            other.b,
            other.c
        ];
        const distsToThisPlane = [
            d1,
            d2,
            d3
        ];
        const ptsOnOther = other._linePlaneIntersect(this.normal, this.d, triVerts2, distsToThisPlane, eps);
        // If we didn’t get at least two intersection points on either triangle, there’s no proper intersection segment
        if (ptsOnThis.length < 2 || ptsOnOther.length < 2) {
            return [];
        }
        // 4) Both sets of intersection points (ptsOnThis and ptsOnOther) lie on the same line (the line of intersection of the two planes).
        // We must compute the overlap of these two segments on that line.
        // Pick two distinct points from ptsOnThis to define the direction of the intersection line
        const p11 = ptsOnThis[0];
        const p12 = ptsOnThis[1];
        // Compute a unit direction vector along the intersection line:
        // direction = normalize( this.normal × other.normal )
        let lineDir = Triangle.cross(this.normal, other.normal);
        lineDir = Triangle.normalize(lineDir);
        // Use p11 as an “origin” on that line. Then we will project all intersection points onto lineDir.
        const origin = p11;
        // Project p11 and p12 onto the line, to get this triangle’s parametric interval [t1Min, t1Max]
        const t11 = Triangle.dot(Triangle.subtract(p11, origin), lineDir);
        const t12 = Triangle.dot(Triangle.subtract(p12, origin), lineDir);
        let t1Min = Math.min(t11, t12);
        let t1Max = Math.max(t11, t12);
        // Similarly project other’s intersection points
        const p21 = ptsOnOther[0];
        const p22 = ptsOnOther[1];
        const t21 = Triangle.dot(Triangle.subtract(p21, origin), lineDir);
        const t22 = Triangle.dot(Triangle.subtract(p22, origin), lineDir);
        let t2Min = Math.min(t21, t22);
        let t2Max = Math.max(t21, t22);
        // The overlap on the line is [max(t1Min, t2Min), min(t1Max, t2Max)]
        const iMin = Math.max(t1Min, t2Min);
        const iMax = Math.min(t1Max, t2Max);
        // If the intervals do not overlap (with a small tolerance), no intersection
        if (iMax < iMin - eps) {
            return [];
        }
        // If they overlap at a single point (within eps), return that single point
        if (Math.abs(iMax - iMin) < eps) {
            const singlePt = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMin));
            return [singlePt];
        }
        // Otherwise they overlap in a segment: return the two endpoints
        const intersectPt1 = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMin));
        const intersectPt2 = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMax));
        return [
            intersectPt1,
            intersectPt2
        ];
    }
    /**
   * Given a triangle’s vertices and their signed distances to a plane (planeNormal · X + planeD = 0),
   * returns the 2 or 3 intersection points between the triangle’s edges and that plane.
   * @param {{x: number, y: number, z: number}} planeNormal
   * @param {number} planeD
   * @param {Array<{x: number,y: number,z: number}>} verts - [v0, v1, v2]
   * @param {Array<number>} dists - [d0, d1, d2] signed distances of each vert to the plane
   * @param {number} eps
   * @returns {Array<{x: number,y: number,z: number}>} 0, 1, or 2 (or 3 if a vertex is exactly on the plane) points
   */
    /**
   * Given a triangle’s vertices and their signed distances to a plane (planeNormal · X + planeD = 0),
   * returns the 2 or 3 intersection points between the triangle’s edges and that plane.
   * @param {{x:number,y:number,z:number}} planeNormal
   * @param {number} planeD
   * @param {Array<{x:number,y:number,z:number}>} verts - [v0, v1, v2]
   * @param {Array<number>} dists - [d0, d1, d2] signed distances of each vert to the plane
   * @param {number} eps
   * @returns {Array<{x:number,y:number,z:number}>} 0, 1, or 2 (or 3 if a vertex is exactly on the plane) points
   */
    _linePlaneIntersect(planeNormal, planeD, verts, dists, eps) {
        const points = [];
        // For each edge (i, j = (i+1)%3):
        for (let i = 0; i < 3; i++) {
            const j = (i + 1) % 3;
            const v1 = verts[i];
            const v2 = verts[j];
            const dist1 = dists[i];
            const dist2 = dists[j];
            // If v1 is very close to the plane, include it
            if (Math.abs(dist1) < eps) {
                points.push({
                    x: v1.x,
                    y: v1.y,
                    z: v1.z
                });
            }
            // If the edge crosses the plane (distances have opposite signs), compute intersection
            if (dist1 * dist2 < -eps) {
                const t = dist1 / (dist1 - dist2);
                const dir = Triangle.subtract(v2, v1);
                const intersectionPoint = Triangle.add(v1, Triangle.multiplyScalar(dir, t));
                points.push(intersectionPoint);
            }
        }
        // Remove duplicates (within eps)
        return Triangle._uniquePoints(points, eps);
    }
    /**
   * Handles the special case where two triangles are (nearly) coplanar.
   * Projects both triangles to 2D (dropping the largest‐magnitude normal component),
   * computes the intersection polygon (convex) via standard segment‐clipping,
   * then lifts the intersection points back to 3D.
   *
   * @param {Triangle} other
   * @param {number} eps
   * @returns {Array<{x: number, y: number, z: number}>} 0 or more points of the overlapped polygon
   */
    /**
   * Handles the special case where two triangles are (nearly) coplanar.
   * Projects both triangles to 2D (dropping the largest‐magnitude normal component),
   * computes the intersection polygon (convex) via standard segment‐clipping,
   * then lifts the intersection points back to 3D.
   *
   * @param {Triangle} other
   * @param {number} eps
   * @returns {Array<{x:number, y:number, z:number}>} 0 or more points of the overlapped polygon
   */
    _coplanarIntersection(other, eps) {
        // 1) Choose which coordinate to “drop” in 3D so we can work in 2D.
        const n = this.normal;
        const ax = Math.abs(n.x), ay = Math.abs(n.y), az = Math.abs(n.z);
        let dropCoord;
        if (ax > ay && ax > az) {
            dropCoord = 'x';
        } else if (ay > az) {
            dropCoord = 'y';
        } else {
            dropCoord = 'z';
        }
        // 2) Project a 3D point to 2D by dropping dropCoord:
        function project2D(pt) {
            if (dropCoord === 'x') {
                return {
                    u: pt.y,
                    v: pt.z
                };
            } else if (dropCoord === 'y') {
                return {
                    u: pt.x,
                    v: pt.z
                };
            } else {
                // dropCoord === 'z'
                return {
                    u: pt.x,
                    v: pt.y
                };
            }
        }
        // Build 2D representations of each triangle
        const tri1_2D = [
            project2D(this.a),
            project2D(this.b),
            project2D(this.c)
        ];
        const tri2_2D = [
            project2D(other.a),
            project2D(other.b),
            project2D(other.c)
        ];
        // Helper: sign of the 2D area (for barycentric test)
        function sign2D(p1, p2, p3) {
            return (p1.u - p3.u) * (p2.v - p3.v) - (p2.u - p3.u) * (p1.v - p3.v);
        }
        // Helper: test if 2D point p is inside the 2D triangle tri (array of 3 points)
        function pointInTri2D(p, tri) {
            const d1 = sign2D(p, tri[0], tri[1]);
            const d2 = sign2D(p, tri[1], tri[2]);
            const d3 = sign2D(p, tri[2], tri[0]);
            const hasNeg = d1 < -eps || d2 < -eps || d3 < -eps;
            const hasPos = d1 > eps || d2 > eps || d3 > eps;
            return !(hasNeg && hasPos);
        }
        // 3) Collect all intersection candidates in 2D:
        const points2D = [];
        // 3a) Any vertex of tri1 that lies inside tri2
        for (let p of tri1_2D) {
            if (pointInTri2D(p, tri2_2D)) {
                points2D.push(p);
            }
        }
        // 3b) Any vertex of tri2 that lies inside tri1
        for (let p of tri2_2D) {
            if (pointInTri2D(p, tri1_2D)) {
                points2D.push(p);
            }
        }
        // 3c) Any intersection between edges of tri1 and edges of tri2
        function segmentIntersection2D(p1, p2, p3, p4) {
            // Solve for intersection of segment p1→p2 and p3→p4 in 2D
            const s1x = p2.u - p1.u;
            const s1y = p2.v - p1.v;
            const s2x = p4.u - p3.u;
            const s2y = p4.v - p3.v;
            const denom = -s2x * s1y + s1x * s2y;
            if (Math.abs(denom) < eps)
                return null;
            // parallel or nearly parallel
            const s = (-s1y * (p1.u - p3.u) + s1x * (p1.v - p3.v)) / denom;
            const t = (s2x * (p1.v - p3.v) - s2y * (p1.u - p3.u)) / denom;
            if (s >= -eps && s <= 1 + eps && t >= -eps && t <= 1 + eps) {
                // They intersect; compute intersection point
                return {
                    u: p1.u + t * s1x,
                    v: p1.v + t * s1y
                };
            }
            return null;
        }
        for (let i = 0; i < 3; i++) {
            const iNext = (i + 1) % 3;
            for (let j = 0; j < 3; j++) {
                const jNext = (j + 1) % 3;
                const inter = segmentIntersection2D(tri1_2D[i], tri1_2D[iNext], tri2_2D[j], tri2_2D[jNext]);
                if (inter) {
                    points2D.push(inter);
                }
            }
        }
        // 4) Deduplicate 2D points (within eps)
        const unique2D = [];
        for (let p of points2D) {
            const found = unique2D.some(q => Math.hypot(p.u - q.u, p.v - q.v) < eps);
            if (!found) {
                unique2D.push(p);
            }
        }
        // 5) “Lift” each 2D point back to 3D by solving for the dropped coordinate from the plane equation:
        //    n.x * x + n.y * y + n.z * z + d = 0
        // If dropCoord = 'x', then x = -(n.y*y + n.z*z + d)/n.x; etc.
        const result3D = [];
        for (let p of unique2D) {
            let x, y, z;
            if (dropCoord === 'x') {
                y = p.u;
                z = p.v;
                // n.x * x + n.y * y + n.z * z + d = 0  =>  x = -(n.y*y + n.z*z + d)/n.x
                x = -(n.y * y + n.z * z + this.d) / n.x;
            } else if (dropCoord === 'y') {
                x = p.u;
                z = p.v;
                y = -(n.x * x + n.z * z + this.d) / n.y;
            } else {
                // dropCoord === 'z'
                x = p.u;
                y = p.v;
                z = -(n.x * x + n.y * y + this.d) / n.z;
            }
            result3D.push({
                x,
                y,
                z
            });
        }
        return result3D;
    }
    // ----------------------- Static Helper Methods -----------------------
    /** Subtract two vectors: a - b */
    static subtract(a, b) {
        return {
            x: a.x - b.x,
            y: a.y - b.y,
            z: a.z - b.z
        };
    }
    /** Add two vectors: a + b */
    static add(a, b) {
        return {
            x: a.x + b.x,
            y: a.y + b.y,
            z: a.z + b.z
        };
    }
    /** Dot product: a · b */
    static dot(a, b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    /** Cross product: a × b */
    static cross(a, b) {
        return {
            x: a.y * b.z - a.z * b.y,
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x
        };
    }
    /** Multiply vector v by scalar s */
    static multiplyScalar(v, s) {
        return {
            x: v.x * s,
            y: v.y * s,
            z: v.z * s
        };
    }
    /** Vector length (magnitude) */
    static length(v) {
        return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }
    /** Normalize vector v to unit length (if nonzero) */
    static normalize(v) {
        const len = Triangle.length(v);
        if (len > 0) {
            return {
                x: v.x / len,
                y: v.y / len,
                z: v.z / len
            };
        }
        return {
            x: 0,
            y: 0,
            z: 0
        };
    }
    /**
   * Remove duplicate 3D points (within epsilon distance).
   * @param {Array<{x: number,y: number,z: number}>} points
   * @param {number} eps
   * @returns {Array<{x: number,y: number,z: number}>}
   */
    /**
   * Remove duplicate 3D points (within epsilon distance).
   * @param {Array<{x:number,y:number,z:number}>} points
   * @param {number} eps
   * @returns {Array<{x:number,y:number,z:number}>}
   */
    static _uniquePoints(points, eps) {
        const unique = [];
        for (let p of points) {
            const found = unique.some(q => Math.hypot(p.x - q.x, p.y - q.y, p.z - q.z) < eps);
            if (!found) {
                unique.push(p);
            }
        }
        return unique;
    }
    /**
   * Split this triangle by the intersection with another triangle.
   * If the two triangles intersect in exactly two points (a line segment),
   * this method subdivides this triangle into smaller triangles along that segment.
   * For any other intersection case (0, 1, or ≥3 points), it returns [this] unchanged.
   *
   * @param {Triangle} other - the other triangle to intersect with
   * @param {number} [eps=1e-6] - epsilon tolerance
   * @returns {Array<Triangle>} An array of new Triangle objects (subdivisions), or [this] if no split is needed.
   */
    splitByTriangle(other, eps = 0.000001) {
        // Compute intersection points
        const pts = this.intersect(other, eps);
        // We only handle the case where there are exactly 2 intersection points (non‐coplanar line intersection).
        if (pts.length !== 2) {
            // No proper line to split, return the original triangle
            return [this];
        }
        // Let P and Q be the two intersection points
        const [P, Q] = pts;
        // Subdivide this triangle by the segment P–Q
        return this._splitByTwoPoints(P, Q, eps);
    }
    /**
   * Internal helper: split *this* triangle by a line segment joining P and Q.
   * Assumes P and Q each lie on one of the edges (or at a vertex) of this triangle, and that they are distinct.
   *
   * @param {{x:number, y:number, z:number}} P - first point on the boundary
   * @param {{x:number, y:number, z:number}} Q - second point on the boundary
   * @param {number} eps - tolerance
   * @returns {Array<Triangle>} Sub‐triangles that exactly cover the original triangle split by P–Q
   */
    _splitByTwoPoints(P, Q, eps) {
        // Original vertices
        const V = [
            this.a,
            this.b,
            this.c
        ];
        // 1) Determine which edge each point lies on.
        const eP = this._findEdgeIndex(P, eps);
        const eQ = this._findEdgeIndex(Q, eps);
        // If either edge index is invalid, just return the original triangle
        if (eP < 0 || eQ < 0) {
            return [this];
        }
        // If both points are on the same edge (or coincide), splitting does nothing.
        if (eP === eQ) {
            return [this];
        }
        // 2) Build two boundary chains (polygons) from P to Q around the triangle edges:
        const poly1 = this._getBoundaryChain(eP, P, eQ, Q, eps);
        const poly2 = this._getBoundaryChainReverse(eP, P, eQ, Q, eps);
        // 3) Triangulate each polygon by "fan" triangulation from the first vertex
        // Collect resulting Triangle objects
        const newTriangles = [];
        // Helper to forge Triangle instances
        const makeTriangle = (p0, p1, p2) => new Triangle(p0, p1, p2);
        // (a) Polygon 1
        for (let k = 1; k + 1 < poly1.length; k++) {
            const p0 = poly1[0];
            const p1 = poly1[k];
            const p2 = poly1[k + 1];
            if (!Triangle._isDegenerate(p0, p1, p2, eps)) {
                newTriangles.push(makeTriangle(p0, p1, p2));
            }
        }
        // (b) Polygon 2
        for (let k = 1; k + 1 < poly2.length; k++) {
            const p0 = poly2[0];
            const p1 = poly2[k];
            const p2 = poly2[k + 1];
            if (!Triangle._isDegenerate(p0, p1, p2, eps)) {
                newTriangles.push(makeTriangle(p0, p1, p2));
            }
        }
        // If no valid subdivisions were produced (e.g. P=Q or degeneracy), return original
        if (newTriangles.length === 0) {
            return [this];
        }
        return newTriangles;
    }
    /**
   * Determine which edge of this triangle a given point lies on.
   * Edges are indexed 0,1,2 corresponding to [a→b], [b→c], [c→a].
   * Returns -1 if the point is not on any edge (within eps).
   *
   * @param {{x:number,y:number,z:number}} pt - the point to test
   * @param {number} eps
   * @returns {number} edge index (0,1,2) or -1 if none
   */
    _findEdgeIndex(pt, eps) {
        const verts = [
            this.a,
            this.b,
            this.c
        ];
        for (let i = 0; i < 3; i++) {
            const v1 = verts[i];
            const v2 = verts[(i + 1) % 3];
            // Check if pt is within eps of the segment v1→v2
            // 1) Compute projection parameter t = ((pt - v1) · (v2 - v1)) / |v2 - v1|^2
            const v1ToV2 = Triangle.subtract(v2, v1);
            const v1ToPt = Triangle.subtract(pt, v1);
            const denom = Triangle.dot(v1ToV2, v1ToV2);
            if (denom < eps)
                continue;
            // edge degeneracy
            const t = Triangle.dot(v1ToPt, v1ToV2) / denom;
            if (t < -eps || t > 1 + eps)
                continue;
            // not between v1 and v2
            // Closest point on infinite line is v1 + t * (v1ToV2)
            const closest = Triangle.add(v1, Triangle.multiplyScalar(v1ToV2, t));
            // Check distance from pt to that closest point
            const dist = Triangle.length(Triangle.subtract(pt, closest));
            if (dist < eps) {
                return i;
            }
        }
        return -1;
    }
    /**
   * Build a boundary chain (one of the two arcs) from P to Q going "forward" around the triangle edges.
   *
   * @param {number} eStart - edge index where P lies
   * @param {{x,y,z}} P - starting point
   * @param {number} eEnd - edge index where Q lies
   * @param {{x,y,z}} Q - ending point
   * @param {number} eps
   * @returns {Array<{x,y,z}>} Array of points [P, (some original vertices), Q]
   */
    _getBoundaryChain(eStart, P, eEnd, Q, eps) {
        const V = [
            this.a,
            this.b,
            this.c
        ];
        const chain = [];
        chain.push(P);
        let currentEdge = eStart;
        let currentVertexIndex = (eStart + 1) % 3;
        // Traverse forward until we reach eEnd
        while (currentEdge !== eEnd) {
            const Vcurr = V[currentVertexIndex];
            // Skip if Vcurr coincides with P or Q (within eps)
            if (Triangle.length(Triangle.subtract(Vcurr, P)) > eps && Triangle.length(Triangle.subtract(Vcurr, Q)) > eps) {
                chain.push(Vcurr);
            }
            currentEdge = (currentEdge + 1) % 3;
            currentVertexIndex = (currentVertexIndex + 1) % 3;
        }
        chain.push(Q);
        return chain;
    }
    /**
   * Build the other boundary chain from P to Q going "backward" around the triangle edges.
   *
   * @param {number} eStart - edge index where P lies
   * @param {{x,y,z}} P - starting point
   * @param {number} eEnd - edge index where Q lies
   * @param {{x,y,z}} Q - ending point
   * @param {number} eps
   * @returns {Array<{x,y,z}>} Array of points [P, (some original vertices), Q]
   */
    _getBoundaryChainReverse(eStart, P, eEnd, Q, eps) {
        const V = [
            this.a,
            this.b,
            this.c
        ];
        const chain = [];
        chain.push(P);
        let currentEdge = eStart;
        let currentVertexIndex = eStart;
        // the "start" vertex for edge eStart is v[eStart]
        // Traverse backward until we reach eEnd
        while (currentEdge !== eEnd) {
            const Vcurr = V[currentVertexIndex];
            // Skip if Vcurr coincides with P or Q
            if (Triangle.length(Triangle.subtract(Vcurr, P)) > eps && Triangle.length(Triangle.subtract(Vcurr, Q)) > eps) {
                chain.push(Vcurr);
            }
            // Step backward: edge - 1 mod 3
            currentEdge = (currentEdge + 2) % 3;
            currentVertexIndex = (currentVertexIndex + 2) % 3;
        }
        chain.push(Q);
        return chain;
    }
    /**
   * Determine if three points form a degenerate (zero‐area) triangle.
   * @param {{x,y,z}} p0
   * @param {{x,y,z}} p1
   * @param {{x,y,z}} p2
   * @param {number} eps
   * @returns {boolean} True if area is < eps
   */
    static _isDegenerate(p0, p1, p2, eps) {
        const v0 = Triangle.subtract(p1, p0);
        const v1 = Triangle.subtract(p2, p0);
        const cross = Triangle.cross(v0, v1);
        const area2 = Triangle.length(cross);
        return area2 < eps;
    }
}