
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
     * @param {{x:number, y:number, z:number}} a - first vertex
     * @param {{x:number, y:number, z:number}} b - second vertex
     * @param {{x:number, y:number, z:number}} c - third vertex
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
        // Unnormalized normal vector
        this.normal = Triangle.cross(Triangle.subtract(b, a), Triangle.subtract(c, a));
        // Plane constant d for plane equation n·X + d = 0
        this.d = -Triangle.dot(this.normal, a);
    }
    /**
     * Compute the intersection between this triangle and another triangle.
     * @param {Triangle} other - the other triangle
     * @param {number} [eps=1e-6] - epsilon tolerance for numerical comparisons
     * @returns {Array<{x:number, y:number, z:number}>} An array of 0, 1, 2, or more points.
     */
    /**
   * Compute intersection points between this triangle and another.
   * Returns [] | [point] | [point, point] | [polygon vertices…]
   */
    intersect(other, eps = 0.000001) {
        // distances of other’s verts to this plane
        const d1 = Triangle.dot(this.normal, other.a) + this.d;
        const d2 = Triangle.dot(this.normal, other.b) + this.d;
        const d3 = Triangle.dot(this.normal, other.c) + this.d;
        if (d1 > eps && d2 > eps && d3 > eps || d1 < -eps && d2 < -eps && d3 < -eps) {
            return [];
        }
        // Coplanar?
        if (Math.abs(d1) < eps && Math.abs(d2) < eps && Math.abs(d3) < eps) {
            return this._coplanarIntersection(other, eps);
        }
        // distances of this verts to other’s plane
        const e1 = Triangle.dot(other.normal, this.a) + other.d;
        const e2 = Triangle.dot(other.normal, this.b) + other.d;
        const e3 = Triangle.dot(other.normal, this.c) + other.d;
        if (e1 > eps && e2 > eps && e3 > eps || e1 < -eps && e2 < -eps && e3 < -eps) {
            return [];
        }
        // find the two points where this triangle cuts other’s plane
        const triVerts1 = [
            this.a,
            this.b,
            this.c
        ];
        const ptsOnThis = this._linePlaneIntersect(other.normal, other.d, triVerts1, [
            e1,
            e2,
            e3
        ], eps);
        if (ptsOnThis.length < 2)
            return [];
        // same for other triangle into this plane
        const triVerts2 = [
            other.a,
            other.b,
            other.c
        ];
        const ptsOnOther = other._linePlaneIntersect(this.normal, this.d, triVerts2, [
            d1,
            d2,
            d3
        ], eps);
        if (ptsOnOther.length < 2)
            return [];
        // Overlap the two segments
        const [p11, p12] = ptsOnThis;
        let lineDir = Triangle.normalize(Triangle.cross(this.normal, other.normal));
        const origin = p11;
        const projectT = p => Triangle.dot(Triangle.subtract(p, origin), lineDir);
        let [t11, t12] = [
            projectT(p11),
            projectT(p12)
        ].sort((a, b) => a - b);
        let [t21, t22] = [
            projectT(ptsOnOther[0]),
            projectT(ptsOnOther[1])
        ].sort((a, b) => a - b);
        const iMin = Math.max(t11, t21), iMax = Math.min(t12, t22);
        if (iMax < iMin - eps)
            return [];
        if (Math.abs(iMax - iMin) < eps) {
            const single = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMin));
            return [single];
        }
        const ip1 = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMin));
        const ip2 = Triangle.add(origin, Triangle.multiplyScalar(lineDir, iMax));
        return [
            ip1,
            ip2
        ];
    }
    /**
     * Given a triangle’s vertices and their signed distances to a plane,
     * returns the intersection points between the triangle’s edges and that plane.
     */
    _linePlaneIntersect(planeNormal, planeD, verts, dists, eps) {
        const points = [];
        for (let i = 0; i < 3; i++) {
            const j = (i + 1) % 3;
            const v1 = verts[i], v2 = verts[j];
            const dist1 = dists[i], dist2 = dists[j];
            if (Math.abs(dist1) < eps) {
                points.push({
                    x: v1.x,
                    y: v1.y,
                    z: v1.z
                });
            }
            if (dist1 * dist2 < -eps) {
                const t = dist1 / (dist1 - dist2);
                const dir = Triangle.subtract(v2, v1);
                const ip = Triangle.add(v1, Triangle.multiplyScalar(dir, t));
                points.push(ip);
            }
        }
        return Triangle._uniquePoints(points, eps);
    }
    /**
     * Handles the special case where two triangles are (nearly) coplanar.
     * Projects to 2D, computes polygon intersection, then lifts back to 3D.
     */
    _coplanarIntersection(other, eps) {
        const n = this.normal;
        const ax = Math.abs(n.x), ay = Math.abs(n.y), az = Math.abs(n.z);
        let drop = ax > ay && ax > az ? 'x' : ay > az ? 'y' : 'z';
        const tri1 = [
            project2D(this.a),
            project2D(this.b),
            project2D(this.c)
        ];
        const tri2 = [
            project2D(other.a),
            project2D(other.b),
            project2D(other.c)
        ];
        function inside(p, tri) {
            const d1 = sign2D(p, tri[0], tri[1]);
            const d2 = sign2D(p, tri[1], tri[2]);
            const d3 = sign2D(p, tri[2], tri[0]);
            const neg = d1 < -eps || d2 < -eps || d3 < -eps;
            const pos = d1 > eps || d2 > eps || d3 > eps;
            return !(neg && pos);
        }
        function segInter(p1, p2, p3, p4) {
            const s1x = p2.u - p1.u, s1y = p2.v - p1.v;
            const s2x = p4.u - p3.u, s2y = p4.v - p3.v;
            const denom = -s2x * s1y + s1x * s2y;
            if (Math.abs(denom) < eps)
                return null;
            const s = (-s1y * (p1.u - p3.u) + s1x * (p1.v - p3.v)) / denom;
            const t = (s2x * (p1.v - p3.v) - s2y * (p1.u - p3.u)) / denom;
            if (s >= -eps && s <= 1 + eps && t >= -eps && t <= 1 + eps) {
                return {
                    u: p1.u + t * s1x,
                    v: p1.v + t * s1y
                };
            }
            return null;
        }
        let points2D = [];
        tri1.forEach(p => inside(p, tri2) && points2D.push(p));
        tri2.forEach(p => inside(p, tri1) && points2D.push(p));
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                const inter = segInter(tri1[i], tri1[(i + 1) % 3], tri2[j], tri2[(j + 1) % 3]);
                if (inter)
                    points2D.push(inter);
            }
        }
        // Deduplicate
        const unique2D = [];
        points2D.forEach(p => {
            if (!unique2D.some(q => Math.hypot(p.u - q.u, p.v - q.v) < eps)) {
                unique2D.push(p);
            }
        });
        // Lift back to 3D
        const result3D = [];
        unique2D.forEach(p => {
            let x, y, z;
            if (drop === 'x') {
                y = p.u;
                z = p.v;
                x = -(n.y * y + n.z * z + this.d) / n.x;
            } else if (drop === 'y') {
                x = p.u;
                z = p.v;
                y = -(n.x * x + n.z * z + this.d) / n.y;
            } else {
                x = p.u;
                y = p.v;
                z = -(n.x * x + n.y * y + this.d) / n.z;
            }
            result3D.push({
                x,
                y,
                z
            });
        });
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
        return len > 0 ? {
            x: v.x / len,
            y: v.y / len,
            z: v.z / len
        } : {
            x: 0,
            y: 0,
            z: 0
        };
    }
    /**
     * Remove duplicate 3D points (within epsilon distance).
     * @param {Array<{x:number,y:number,z:number}>} points
     * @param {number} eps
     * @returns {Array<{x:number,y:number,z:number}>}
     */
    static _uniquePoints(points, eps) {
        const unique = [];
        points.forEach(p => {
            if (!unique.some(q => Math.hypot(p.x - q.x, p.y - q.y, p.z - q.z) < eps)) {
                unique.push(p);
            }
        });
        return unique;
    }
    /**
   * Split this triangle by the intersection with otherTriangle.
   * If there’s a proper line-segment intersection, returns an array of
   * smaller Triangle instances whose edges include that segment.
   * If no segment (or coplanar), returns [this].
   */
    splitByTriangle(other, eps = 0.000001) {
        const d_other_a_on_this_plane = Triangle.dot(this.normal, other.a) + this.d;
        const d_other_b_on_this_plane = Triangle.dot(this.normal, other.b) + this.d;
        const d_other_c_on_this_plane = Triangle.dot(this.normal, other.c) + this.d;
        const is_coplanar = Math.abs(d_other_a_on_this_plane) < eps && Math.abs(d_other_b_on_this_plane) < eps && Math.abs(d_other_c_on_this_plane) < eps;
        if (is_coplanar) {
            // Coplanar case (unchanged from previous correct version)
            let current_triangles = [this];
            const other_edges = [
                {
                    v1: other.a,
                    v2: other.b
                },
                {
                    v1: other.b,
                    v2: other.c
                },
                {
                    v1: other.c,
                    v2: other.a
                }
            ];
            for (const edge of other_edges) {
                const edge_vec = Triangle.subtract(edge.v2, edge.v1);
                if (Triangle.length(edge_vec) < eps)
                    continue;
                const line_normal_3d = Triangle.normalize(Triangle.cross(edge_vec, this.normal));
                if (Triangle.length(line_normal_3d) < eps || Triangle.length(this.normal) < eps) {
                    continue;
                }
                const line_d_param = -Triangle.dot(line_normal_3d, edge.v1);
                let next_triangles_batch = [];
                for (const current_tri of current_triangles) {
                    const tv = [
                        current_tri.a,
                        current_tri.b,
                        current_tri.c
                    ];
                    const dists = tv.map(v_node => Triangle.dot(line_normal_3d, v_node) + line_d_param);
                    const all_positive_side = dists.every(d_node => d_node > eps);
                    const all_negative_side = dists.every(d_node => d_node < -eps);
                    const all_on_or_near_line = dists.every(d_node => Math.abs(d_node) < eps);
                    if (all_positive_side || all_negative_side || all_on_or_near_line) {
                        next_triangles_batch.push(current_tri);
                        continue;
                    }
                    const cut_points = current_tri._linePlaneIntersect(line_normal_3d, line_d_param, tv, dists, eps);
                    if (cut_points.length === 2) {
                        const [p_cut, q_cut] = cut_points;
                        if (Triangle.length(Triangle.subtract(p_cut, q_cut)) > eps) {
                            next_triangles_batch.push(...current_tri._splitBySegment(p_cut, q_cut, eps));
                        } else {
                            next_triangles_batch.push(current_tri);
                        }
                    } else {
                        next_triangles_batch.push(current_tri);
                    }
                }
                current_triangles = next_triangles_batch;
            }
            return current_triangles.filter(t => {
                // Filter degenerate final triangles
                const edge1 = Triangle.subtract(t.b, t.a);
                const edge2 = Triangle.subtract(t.c, t.a);
                return Triangle.length(Triangle.cross(edge1, edge2)) > eps * eps;
            });
        } else {
            // Non-coplanar case
            const actual_intersection_points = this.intersect(other, eps);
            if (actual_intersection_points.length !== 2) {
                return [this];
            }
            // No segment intersection
            let [p_int, q_int] = actual_intersection_points;
            const point_equals = (pt1, pt2, eq_eps) => Math.hypot(pt1.x - pt2.x, pt1.y - pt2.y, pt1.z - pt2.z) < eq_eps;
            const is_point_on_segment = (pt, seg_a, seg_b, on_eps) => {
                if (point_equals(pt, seg_a, on_eps) || point_equals(pt, seg_b, on_eps))
                    return true;
                // Point is one of the segment endpoints
                const vec_sa_sb = Triangle.subtract(seg_b, seg_a);
                const vec_sa_pt = Triangle.subtract(pt, seg_a);
                const len_sq_sa_sb = Triangle.dot(vec_sa_sb, vec_sa_sb);
                if (len_sq_sa_sb < on_eps * on_eps)
                    return false;
                // Segment is degenerate
                const cross_prod = Triangle.cross(vec_sa_pt, vec_sa_sb);
                const dist_sq_pt_line = Triangle.dot(cross_prod, cross_prod) / len_sq_sa_sb;
                if (dist_sq_pt_line > on_eps * on_eps)
                    return false;
                // Point not collinear
                const dot_prod = Triangle.dot(vec_sa_pt, vec_sa_sb);
                // Check if projection is within segment [seg_a, seg_b]
                return dot_prod >= -on_eps * on_eps && dot_prod <= len_sq_sa_sb + on_eps * on_eps;
            };
            const get_point_status = (pt, tri_to_check, status_eps) => {
                const vertices = [
                    tri_to_check.a,
                    tri_to_check.b,
                    tri_to_check.c
                ];
                for (let i = 0; i < 3; ++i) {
                    if (point_equals(pt, vertices[i], status_eps))
                        return 'vertex';
                }
                for (let i = 0; i < 3; ++i) {
                    if (is_point_on_segment(pt, vertices[i], vertices[(i + 1) % 3], status_eps))
                        return 'edge';
                }
                // Check if pt is on plane of tri_to_check (scaled tolerance)
                const tri_normal_len = Triangle.length(tri_to_check.normal);
                if (tri_normal_len < status_eps)
                    return 'degenerate_triangle_plane';
                // Cannot determine for degenerate triangle
                if (Math.abs(Triangle.dot(tri_to_check.normal, pt) + tri_to_check.d) > status_eps * tri_normal_len) {
                    return 'outside_plane';
                }
                // Should not happen if pt from this.intersect()
                // Project to 2D and use barycentric/sign test (like in _coplanarIntersection's 'inside')
                const ax = Math.abs(tri_to_check.normal.x), ay = Math.abs(tri_to_check.normal.y), az = Math.abs(tri_to_check.normal.z);
                let drop = ax > ay && ax > az ? 'x' : ay > az ? 'y' : 'z';
                function project2D(p3d) {
                    if (drop === 'x')
                        return {
                            u: p3d.y,
                            v: p3d.z
                        };
                    if (drop === 'y')
                        return {
                            u: p3d.x,
                            v: p3d.z
                        };
                    return {
                        u: p3d.x,
                        v: p3d.y
                    };
                }
                const p2d = project2D(pt);
                const tri2d_verts = [
                    project2D(tri_to_check.a),
                    project2D(tri_to_check.b),
                    project2D(tri_to_check.c)
                ];
                function sign2D(p1, p2, p3) {
                    return (p1.u - p3.u) * (p2.v - p3.v) - (p2.u - p3.u) * (p1.v - p3.v);
                }
                const s1 = sign2D(p2d, tri2d_verts[0], tri2d_verts[1]);
                const s2 = sign2D(p2d, tri2d_verts[1], tri2d_verts[2]);
                const s3 = sign2D(p2d, tri2d_verts[2], tri2d_verts[0]);
                const has_neg = s1 < -status_eps || s2 < -status_eps || s3 < -status_eps;
                const has_pos = s1 > status_eps || s2 > status_eps || s3 > status_eps;
                if (!(has_neg && has_pos))
                    return 'interior';
                // On or inside projected triangle
                return 'outside';
            };
            // On plane, but outside triangle bounds
            let status_p = get_point_status(p_int, this, eps);
            let status_q = get_point_status(q_int, this, eps);
            const V = [
                this.a,
                this.b,
                this.c
            ];
            let result_triangles = [];
            // Normalize statuses: if outside_plane or degenerate_triangle_plane, treat as error/no_split
            if (status_p === 'outside_plane' || status_q === 'outside_plane' || status_p === 'degenerate_triangle_plane' || status_q === 'degenerate_triangle_plane' || status_p === 'outside' || status_q === 'outside') {
                // Points from intersect should be on/in.
                return [this];
            }
            if ((status_p === 'vertex' || status_p === 'edge') && (status_q === 'vertex' || status_q === 'edge')) {
                // Scenario A: Both points on boundary.
                result_triangles = this._splitBySegment(p_int, q_int, eps);
            } else if (status_p === 'interior' && (status_q === 'vertex' || status_q === 'edge')) {
                // Scenario B: P_int interior, Q_int on boundary.
                for (let i = 0; i < 3; ++i) {
                    if (point_equals(q_int, V[i], eps)) {
                        // Q_int is vertex V[i]
                        result_triangles.push(new Triangle(p_int, V[i], V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(p_int, V[i], V[(i + 2) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 1) % 3], V[(i + 2) % 3]));
                        break;
                    }
                    if (is_point_on_segment(q_int, V[i], V[(i + 1) % 3], eps)) {
                        // Q_int on edge V[i]V[(i+1)%3]
                        result_triangles.push(new Triangle(p_int, V[i], q_int));
                        result_triangles.push(new Triangle(p_int, q_int, V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 1) % 3], V[(i + 2) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 2) % 3], V[i]));
                        break;
                    }
                }
            } else if ((status_p === 'vertex' || status_p === 'edge') && status_q === 'interior') {
                // Scenario B symmetric: Q_int interior, P_int on boundary. Swap them.
                [p_int, q_int] = [
                    q_int,
                    p_int
                ];
                // p_int is now interior, q_int on boundary
                for (let i = 0; i < 3; ++i) {
                    if (point_equals(q_int, V[i], eps)) {
                        // Q_int is vertex V[i]
                        result_triangles.push(new Triangle(p_int, V[i], V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(p_int, V[i], V[(i + 2) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 1) % 3], V[(i + 2) % 3]));
                        break;
                    }
                    if (is_point_on_segment(q_int, V[i], V[(i + 1) % 3], eps)) {
                        // Q_int on edge V[i]V[(i+1)%3]
                        result_triangles.push(new Triangle(p_int, V[i], q_int));
                        result_triangles.push(new Triangle(p_int, q_int, V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 1) % 3], V[(i + 2) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 2) % 3], V[i]));
                        break;
                    }
                }
            } else if (status_p === 'interior' && status_q === 'interior') {
                // Scenario C: Both P_int, Q_int interior.
                let q_found_in_sub_triangle = false;
                for (let i = 0; i < 3; ++i) {
                    const sub_tri_for_q_check = new Triangle(p_int, V[i], V[(i + 1) % 3]);
                    const q_status_in_sub = get_point_status(q_int, sub_tri_for_q_check, eps);
                    if (q_status_in_sub === 'interior' || q_status_in_sub === 'edge' || q_status_in_sub === 'vertex') {
                        result_triangles.push(new Triangle(q_int, p_int, V[i]));
                        result_triangles.push(new Triangle(q_int, p_int, V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(q_int, V[i], V[(i + 1) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 1) % 3], V[(i + 2) % 3]));
                        result_triangles.push(new Triangle(p_int, V[(i + 2) % 3], V[i]));
                        q_found_in_sub_triangle = true;
                        break;
                    }
                }
                if (!q_found_in_sub_triangle) {
                    // Should not happen if P, Q distinct and interior
                    result_triangles = [this];
                }
            } else {
                // Fallback, should not be reached if statuses are comprehensive
                result_triangles = [this];
            }
            if (result_triangles.length === 0)
                return [this];
            // If logic failed to produce any triangles
            return result_triangles.filter(t => {
                const edge1 = Triangle.subtract(t.b, t.a);
                const edge2 = Triangle.subtract(t.c, t.a);
                // Check area; cross product magnitude is 2*Area. Use eps*eps for area-like tolerance.
                return Triangle.length(Triangle.cross(edge1, edge2)) > eps * eps;
            });
        }
    }
    /**
   * Internal: split this triangle by a chord from p→q lying on its edges.
   * Returns an array of Triangle pieces.
   */
    _splitBySegment(p, q, eps) {
        // This method assumes p and q are on the boundary of `this` triangle.
        const same = (u, v) => Math.hypot(u.x - v.x, u.y - v.y, u.z - v.z) < eps;
        const original_boundary_vertices = [
            this.a,
            this.b,
            this.c
        ];
        const insertions = [];
        [
            p,
            q
        ].forEach(pt_to_insert => {
            if (original_boundary_vertices.some(v_orig => same(pt_to_insert, v_orig))) {
                return;
            }
            for (let i = 0; i < 3; i++) {
                const v1 = original_boundary_vertices[i];
                const v2 = original_boundary_vertices[(i + 1) % 3];
                const edge_dir = Triangle.subtract(v2, v1);
                const edge_len_sq = Triangle.dot(edge_dir, edge_dir);
                if (edge_len_sq < eps * eps)
                    continue;
                const t_param = Triangle.dot(Triangle.subtract(pt_to_insert, v1), edge_dir) / edge_len_sq;
                if (t_param > eps && t_param < 1 - eps) {
                    const projected_pt = Triangle.add(v1, Triangle.multiplyScalar(edge_dir, t_param));
                    if (same(pt_to_insert, projected_pt)) {
                        insertions.push({
                            point: pt_to_insert,
                            after_idx: i,
                            t_param: t_param
                        });
                        break;
                    }
                }
            }
        });
        insertions.sort((a, b) => {
            if (a.after_idx !== b.after_idx) {
                return b.after_idx - a.after_idx;
            }
            return b.t_param - a.t_param;
        });
        const augmented_boundary = [...original_boundary_vertices];
        insertions.forEach(({ point, after_idx }) => {
            augmented_boundary.splice(after_idx + 1, 0, point);
        });
        const idxP = augmented_boundary.findIndex(v => same(v, p));
        const idxQ = augmented_boundary.findIndex(v => same(v, q));
        if (idxP === -1 || idxQ === -1 || idxP === idxQ) {
            return [this];
        }
        const n_boundary_pts = augmented_boundary.length;
        const poly1_pts = [], poly2_pts = [];
        for (let i = idxP; ; i = (i + 1) % n_boundary_pts) {
            poly1_pts.push(augmented_boundary[i]);
            if (i === idxQ)
                break;
        }
        for (let i = idxQ; ; i = (i + 1) % n_boundary_pts) {
            poly2_pts.push(augmented_boundary[i]);
            if (i === idxP)
                break;
        }
        const result_triangles = [];
        [
            poly1_pts,
            poly2_pts
        ].forEach(poly_vertices => {
            if (poly_vertices.length < 3)
                return;
            const v0 = poly_vertices[0];
            for (let i = 1; i < poly_vertices.length - 1; i++) {
                const v_i = poly_vertices[i];
                const v_i_plus_1 = poly_vertices[i + 1];
                const edge1 = Triangle.subtract(v_i, v0);
                const edge2 = Triangle.subtract(v_i_plus_1, v0);
                const cross_prod_mag = Triangle.length(Triangle.cross(edge1, edge2));
                if (cross_prod_mag > eps * eps) {
                    // Area-like tolerance
                    result_triangles.push(new Triangle(v0, v_i, v_i_plus_1));
                }
            }
        });
        if (result_triangles.length === 0)
            return [this];
        return result_triangles;
    }
}