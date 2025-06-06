// script.js

import * as THREE from 'three';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls.js';
import { Triangle } from './kernel.js'; // Your Triangle class with splitByTriangle

let scene, camera, renderer, controls;
let tri1Mesh = null,
    tri2Mesh = null;
let tri1Pieces = [],
    tri2Pieces = [];
let intersectGroup = null;

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

    // 7) Hook up the “Update Triangles” button
    document.getElementById('updateBtn').addEventListener('click', () => {
        updateTriangles();
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

// add an onchange event to all input fields to update triangles on change including the split mode checkbox
document.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
    input.addEventListener('change', () => {
        updateTriangles();
    });
});








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

    // Check split mode
    const splitMode = document.getElementById('splitMode').checked;

    if (!splitMode) {
        // DEFAULT: render original triangles
        tri1Mesh = createTriangleMesh([t1a, t1b, t1c], 0x2194ce, 0.5);
        tri2Mesh = createTriangleMesh([t2a, t2b, t2c], 0xce3219, 0.5);
        scene.add(tri1Mesh);
        scene.add(tri2Mesh);
    } else {
        // SPLIT MODE: split each triangle by the other, then render all resulting sub-triangles
        const tri1Subs = tri1.splitByTriangle(tri2, 1e-6);
        const tri2Subs = tri2.splitByTriangle(tri1, 1e-6);


        tri1Subs.forEach(subTri => {
            // generate a random color for each sub-triangle
            const myColor = new THREE.Color(Math.random(), Math.random(), Math.random());

            const pts = [subTri.a, subTri.b, subTri.c];
            const mesh = createTriangleMesh(pts, myColor, 0.5);
            scene.add(mesh);
            tri1Pieces.push(mesh);
        });
        tri2Subs.forEach(subTri => {            
            // generate a random color for each sub-triangle
            const myColor = new THREE.Color(Math.random(), Math.random(), Math.random());

            const pts = [subTri.a, subTri.b, subTri.c];
            const mesh = createTriangleMesh(pts, myColor, 0.5);
            scene.add(mesh);
            tri2Pieces.push(mesh);
        });
    }

    // Compute and render intersection geometry exactly as before
    const intersectionPoints = tri1.intersect(tri2, 1e-6);
    intersectGroup = new THREE.Group();

    if (intersectionPoints.length === 1) {
        const pt = intersectionPoints[0];
        const sphereGeom = new THREE.SphereGeometry(0.05, 16, 16);
        const sphereMat = new THREE.MeshBasicMaterial({ color: 0xffff00 });
        const sphere = new THREE.Mesh(sphereGeom, sphereMat);
        sphere.position.set(pt.x, pt.y, pt.z);
        intersectGroup.add(sphere);

    } else if (intersectionPoints.length === 2) {
        const pts = intersectionPoints.map(p => new THREE.Vector3(p.x, p.y, p.z));
        const lineGeom = new THREE.BufferGeometry().setFromPoints(pts);
        const lineMat = new THREE.LineBasicMaterial({ color: 0xffff00, linewidth: 2 });
        const line = new THREE.Line(lineGeom, lineMat);
        intersectGroup.add(line);

    } else if (intersectionPoints.length > 2) {
        // Coplanar overlap polygon
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
        if (Triangle.length(u) < 1e-6 && intersectionPoints.length > 1) {
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
            const xu = Triangle.dot(rel, u);
            const yv = Triangle.dot(rel, v);
            const angle = Math.atan2(yv, xu);
            return { point: p, angle };
        });
        pointsWithAngle.sort((a, b) => a.angle - b.angle);
        const sortedPoints = pointsWithAngle.map(pa => pa.point);

        const loopVecs = sortedPoints.map(p => new THREE.Vector3(p.x, p.y, p.z));
        loopVecs.push(new THREE.Vector3(sortedPoints[0].x, sortedPoints[0].y, sortedPoints[0].z));
        const loopGeom = new THREE.BufferGeometry().setFromPoints(loopVecs);
        const loopMat = new THREE.LineBasicMaterial({ color: 0xffff00, linewidth: 2 });
        const loop = new THREE.LineLoop(loopGeom, loopMat);
        intersectGroup.add(loop);
    }

    scene.add(intersectGroup);
}

function createTriangleMesh(pts, hexColor, opacity) {
    const geometry = new THREE.BufferGeometry();
    const vertices = new Float32Array([
        pts[0].x, pts[0].y, pts[0].z,
        pts[1].x, pts[1].y, pts[1].z,
        pts[2].x, pts[2].y, pts[2].z
    ]);
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
    geometry.setIndex([0, 1, 2]);
    geometry.computeVertexNormals();

    const material = new THREE.MeshBasicMaterial({
        color: hexColor,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: opacity,
        wireframe: false
    });

    return new THREE.Mesh(geometry, material);
}
