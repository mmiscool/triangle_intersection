// script.js

import * as THREE from 'three';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls.js';
import { Triangle } from './kernel.js'; // Your Triangle class with splitByTriangle and _splitByTwoPoints

let scene, camera, renderer, controls;
let tri1Mesh = null,
    tri2Mesh = null;
let tri1Pieces = [],
    tri2Pieces = [];
let intersectGroup = null;

const EPS = 1e-6;

init();
animate();

function init() {
    // 1) Scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x111111);

    // 2) Orthographic camera
    const aspect = (window.innerWidth - 300) / window.innerHeight;
    const frustumSize = 10;
    camera = new THREE.OrthographicCamera(
        -frustumSize * aspect / 2,
        +frustumSize * aspect / 2,
        +frustumSize / 2,
        -frustumSize / 2,
        0.1,
        1000
    );
    camera.position.set(5, 5, 5);
    camera.lookAt(0, 0, 0);

    // 3) Renderer
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth - 300, window.innerHeight);
    document.getElementById('rendererContainer').appendChild(renderer.domElement);

    // 4) TrackballControls
    controls = new TrackballControls(camera, renderer.domElement);
    controls.rotateSpeed = 4.0;
    controls.zoomSpeed = 2.0;
    controls.panSpeed = 1.0;
    controls.staticMoving = true;

    // 5) Axes helper
    const axesHelper = new THREE.AxesHelper(5);
    scene.add(axesHelper);

    // 6) Hemisphere light (optional)
    const hemiLight = new THREE.HemisphereLight(0xffffff, 0x444444, 1.0);
    hemiLight.position.set(0, 20, 0);
    scene.add(hemiLight);

    // 7) Hook up the "Update Triangles" button and input change
    document.getElementById('updateBtn').addEventListener('click', updateTriangles);
    document.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
        input.addEventListener('change', updateTriangles);
    });

    // 8) Initial draw
    updateTriangles();

    // 9) Window resize handling
    window.addEventListener('resize', onWindowResize);
}

function onWindowResize() {
    const aspect = (window.innerWidth - 300) / window.innerHeight;
    const frustumSize = 10;

    camera.left = -frustumSize * aspect / 2;
    camera.right = +frustumSize * aspect / 2;
    camera.top = +frustumSize / 2;
    camera.bottom = -frustumSize / 2;
    camera.updateProjectionMatrix();

    renderer.setSize(window.innerWidth - 300, window.innerHeight);
    controls.handleResize();
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}

function updateTriangles() {
    // CLEANUP previous meshes/groups
    if (tri1Mesh) {
        scene.remove(tri1Mesh);
        tri1Mesh.geometry.dispose();
        tri1Mesh.material.dispose();
        tri1Mesh = null;
    }
    if (tri2Mesh) {
        scene.remove(tri2Mesh);
        tri2Mesh.geometry.dispose();
        tri2Mesh.material.dispose();
        tri2Mesh = null;
    }
    tri1Pieces.forEach(m => {
        scene.remove(m);
        m.geometry.dispose();
        m.material.dispose();
    });
    tri1Pieces = [];
    tri2Pieces.forEach(m => {
        scene.remove(m);
        m.geometry.dispose();
        m.material.dispose();
    });
    tri2Pieces = [];
    if (intersectGroup) {
        intersectGroup.children.forEach(child => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
        scene.remove(intersectGroup);
        intersectGroup = null;
    }

    // READ input fields for each point
    function readPt(prefix) {
        const x = parseFloat(document.getElementById(prefix + '_x').value);
        const y = parseFloat(document.getElementById(prefix + '_y').value);
        const z = parseFloat(document.getElementById(prefix + '_z').value);
        return { x, y, z };
    }
    const t1a = readPt('t1a');
    const t1b = readPt('t1b');
    const t1c = readPt('t1c');
    const t2a = readPt('t2a');
    const t2b = readPt('t2b');
    const t2c = readPt('t2c');

    // Instantiate Triangle objects
    const tri1 = new Triangle(t1a, t1b, t1c);
    const tri2 = new Triangle(t2a, t2b, t2c);

    // Compute intersection points once
    const intersectionPoints = tri1.intersect(tri2, EPS);

    // Check split mode
    const splitMode = document.getElementById('splitMode').checked;

    if (splitMode && intersectionPoints.length >= 2) {

        // Split triangle 1, assign each piece a unique HSL‐based color
        const pieces1 = tri1.splitByTriangle(tri2, EPS);
        console.log('Split pieces for triangle 1:', pieces1);


        pieces1.forEach((piece, idx) => {
            const h = idx / pieces1.length;                        // vary hue
            const colorHex = new THREE.Color().setHSL(h, 0.7, 0.5).getHex();
            const mesh = createTriangleMesh([piece.a, piece.b, piece.c], colorHex, 0.6);
            tri1Pieces.push(mesh);
            scene.add(mesh);
        });

        // Split triangle 2, continue hue sequence so no two pieces share hue
        const pieces2 = tri2.splitByTriangle(tri1, EPS);
        pieces2.forEach((piece, idx) => {
            const h = (pieces1.length + idx) / (pieces1.length + pieces2.length);
            const colorHex = new THREE.Color().setHSL(h, 0.7, 0.5).getHex();
            const mesh = createTriangleMesh([piece.a, piece.b, piece.c], colorHex, 0.6);
            tri2Pieces.push(mesh);
            scene.add(mesh);
        });




        const debugObject = {
            tri1: tri1,
            tri2: tri2,
            pieces1: pieces1,
            pieces2: pieces2,
            intersectionPoints: intersectionPoints
        }

        console.log('Debug Info:', JSON.stringify(debugObject));

    } else {
        // Default: render original triangles with two distinct colors
        const base1 = new THREE.Color().setHSL(0.0, 0.7, 0.5).getHex();
        const base2 = new THREE.Color().setHSL(0.5, 0.7, 0.5).getHex();

        tri1Mesh = createTriangleMesh([t1a, t1b, t1c], base1, 0.5);
        tri2Mesh = createTriangleMesh([t2a, t2b, t2c], base2, 0.5);
        scene.add(tri1Mesh);
        scene.add(tri2Mesh);
    }

    // Compute and render intersection geometry (unchanged)…
    intersectGroup = new THREE.Group();
    if (intersectionPoints.length === 1) {
        const pt = intersectionPoints[0];
        const sphere = new THREE.Mesh(
            new THREE.SphereGeometry(0.05, 16, 16),
            new THREE.MeshBasicMaterial({ color: 0xffff00 })
        );
        sphere.position.set(pt.x, pt.y, pt.z);
        intersectGroup.add(sphere);
    } else if (intersectionPoints.length === 2) {
        const pts = intersectionPoints.map(p => new THREE.Vector3(p.x, p.y, p.z));
        const line = new THREE.Line(
            new THREE.BufferGeometry().setFromPoints(pts),
            new THREE.LineBasicMaterial({ color: 0xffff00, linewidth: 2 })
        );
        intersectGroup.add(line);
    } else if (intersectionPoints.length > 2) {
        // …coplanar overlap polygon code remains the same…
        const centroid = intersectionPoints.reduce((acc, p) => {
            acc.x += p.x; acc.y += p.y; acc.z += p.z;
            return acc;
        }, { x: 0, y: 0, z: 0 });
        centroid.x /= intersectionPoints.length;
        centroid.y /= intersectionPoints.length;
        centroid.z /= intersectionPoints.length;

        const normal = Triangle.normalize(tri1.normal);
        let v0 = {
            x: intersectionPoints[0].x - centroid.x,
            y: intersectionPoints[0].y - centroid.y,
            z: intersectionPoints[0].z - centroid.z
        };
        let u = Triangle.normalize(v0);
        if (Triangle.length(u) < EPS && intersectionPoints.length > 1) {
            const v1 = {
                x: intersectionPoints[1].x - centroid.x,
                y: intersectionPoints[1].y - centroid.y,
                z: intersectionPoints[1].z - centroid.z
            };
            u = Triangle.normalize(v1);
        }
        const v = Triangle.normalize(Triangle.cross(normal, u));

        const pointsWithAngle = intersectionPoints.map(p => {
            const rel = {
                x: p.x - centroid.x,
                y: p.y - centroid.y,
                z: p.z - centroid.z
            };
            const angle = Math.atan2(
                Triangle.dot(rel, v),
                Triangle.dot(rel, u)
            );
            return { point: p, angle };
        });
        pointsWithAngle.sort((a, b) => a.angle - b.angle);
        const sorted = pointsWithAngle.map(pa => pa.point);

        const loopVecs = sorted.map(p => new THREE.Vector3(p.x, p.y, p.z));
        loopVecs.push(loopVecs[0].clone());
        const loop = new THREE.LineLoop(
            new THREE.BufferGeometry().setFromPoints(loopVecs),
            new THREE.LineBasicMaterial({ color: 0xffff00, linewidth: 2 })
        );
        intersectGroup.add(loop);
    }
    scene.add(intersectGroup);
}


function createTriangleMesh(points, color, opacity) {
    // Build BufferGeometry from 3 points
    const geometry = new THREE.BufferGeometry();
    const positions = new Float32Array(points.length * 3);
    for (let i = 0; i < points.length; i++) {
        positions[3 * i] = points[i].x;
        positions[3 * i + 1] = points[i].y;
        positions[3 * i + 2] = points[i].z;
    }
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setIndex([0, 1, 2]);
    geometry.computeVertexNormals();

    // Create a double-sided, transparent material
    const material = new THREE.MeshBasicMaterial({
        color: color,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: opacity
    });

    // Return the mesh
    return new THREE.Mesh(geometry, material);
}
