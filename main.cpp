/*
 * CSC 305 201801 UVIC
 * The purpose of this source file is to demonstrate the Mesh class which you may use in assignment 2
 * Its only functionality is to render vertices/normals/textures and load textures from png files
 * A demonstration of an ImGui menu window is also included in this file
*/
#include "Mesh/Mesh.h"
#include "OpenGP/GL/glfw_helpers.h"

#include <cstdio>
#include <OpenGP/types.h>
#include <OpenGP/MLogger.h>
#include <OpenGP/GL/Application.h>
#include <OpenGP/GL/ImguiRenderer.h>

using namespace OpenGP;
const float pi = 4*atan(1);

/// Calculates the normal of a triangle in 3-space based on the given 3 points
Vec3 calculateNormal(Vec3 A, Vec3 B){
    return A.cross(B);    //might be backwards; double-check
}

/// Produces a normal for each triangle & adds it to normList
void getNormals(std::vector<Vec3> &vertList, std::vector<unsigned int> & indexList, std::vector<Vec3> &normList){
    for(int i = 0; i < indexList.size(); i+=3){
        Vec3 normal = calculateNormal(vertList[indexList[i]],vertList[indexList[i+1]]);
        normList.push_back(normal);
    }
}

/// Maps the given point to UV-coordinates on a sphere with origin o & radius r
Vec2 sphereMap(Vec3 p, Vec3 o, float r){
    float the, phi, u, v;
    //indices may need to change (x = left->right, y = down->up, z = in->out)
    the = atan2(p[0]-o[0],p[2]-o[2]);
    phi = acos((p[1]-o[1])/r);
    v = phi/2*pi;
    u = (pi-the)/pi;
    return Vec2(-u,-v);
}

/// Loads a cube into renderMesh via vertList, indexList, and normList, with the given origin point & side length
void loadCube(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<Vec2> &tCoordList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, const char * texFile, Vec3 origin, float side){
    std::vector<int> indices = {0,1,2,2,3,0,4,0,3,3,7,4,7,4,5,5,6,7,1,5,6,6,2,1,1,0,4,4,5,1,2,3,7,7,6,2};

    float dx = origin[0]+side;
    float dy = origin[1]+side;
    float dz = origin[2]+side;

    vertList.push_back(origin);                         //0
    vertList.push_back(Vec3(dx,origin[1],origin[2]));   //1
    vertList.push_back(Vec3(dx,dy,origin[2]));          //2
    vertList.push_back(Vec3(origin[0],dy,origin[2]));   //3
    vertList.push_back(Vec3(origin[0],origin[1],dz));   //4
    vertList.push_back(Vec3(dx,origin[1],dz));          //5
    vertList.push_back(Vec3(dx,dy,dz));                 //6
    vertList.push_back(Vec3(origin[0],dy,dz));          //7

    for(int i = 0; i < indices.size(); i++){
        indexList.push_back(indices[i]);
    }

    renderMesh.loadVertices(vertList, indexList);

    getNormals(vertList,indexList,normList);

    renderMesh.loadNormals(normList);

    /// Load textures (assumes texcoords)
    renderMesh.loadTextures(texFile);

    /// Load texture coordinates (assumes textures)

    renderMesh.loadTexCoords(tCoordList);
}

///Returns the midpoint between two points on a triangle
Vec3 getMidpoint(Vec3 A, Vec3 B){
    return Vec3((A[0]+B[0])/2.0,(A[1]+B[1])/2,(A[2]+B[2])/2);   //this doesn't check for duplicates :(
}

/// Loads an icosphere into renderMesh via vertList, indexList, and normList, with the given origin point, radius, & level of subdivision (pending)
/// based off pseudocode from: https://www.csee.umbc.edu/~squire/reference/polyhedra.shtml#icosahedron
/// plus subdivision algorithm from http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
void loadIcoSphere(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<Vec2> &tCoordList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, const char * texFile, Vec3 origin, float radius, int sublevel){
    float vertices[12][3];      //icosahedron x,y,z coordinates
    std::vector<int> indices = {0,1,2,0,2,3,0,3,4,0,4,5,0,5,1,11,6,7,11,7,8,11,8,9,11,9,10,11,10,6,1,2,6,2,3,7,3,4,8,4,5,9,5,1,10,6,7,2,7,8,3,8,9,4,9,10,5,10,6,1};

    /// Angles needed for icosahedron
    float phiaa = 26.56505;
    float phia = pi*phiaa/180.0;
    float theb = pi*36.0/180.0;
    float the72 = pi*72.0/180.0;

    /// Create top & bottom vertices
    vertices[0][0] = 0.0f;
    vertices[0][1] = 0.0f;
    vertices[0][2] = radius;
    vertices[11][0] = 0.0f;
    vertices[11][1] = 0.0f;
    vertices[11][2] = -radius;

    float the = 0.0f;

    /// Generate top half of icosahedron vertices
    for(int i = 1; i < 6; i++){
        vertices[i][0] = radius*cos(the)*cos(phia);
        vertices[i][1] = radius*sin(the)*cos(phia);
        vertices[i][2] = radius*sin(phia);
        the += the72;
    }

    the = theb;
    /// Generate bottom half of icosahedron vertices
    for(int i = 6; i < 11; i++){
        vertices[i][0] = radius*cos(the)*cos(-phia);
        vertices[i][1] = radius*sin(the)*cos(-phia);
        vertices[i][2] = radius*sin(-phia);
        the += the72;
    }

    /// Put vertices into vertList in order
    for(int i = 0; i < 12; i++){
        vertList.push_back(Vec3(vertices[i][0],vertices[i][1],vertices[i][2]));
    }

    /// Put indices into indexList in face order
    for(int i = 0; i < indices.size(); i++){
        indexList.push_back(indices[i]);
    }

    ///TODO: Subdivide each triangle into 4
    //subdivides first triangle, then performs strangely
    for(int i = 0; i < sublevel; i++){
        std::vector<Vec3> tempVertList;
        std::vector<unsigned int> tempIndexList;
        int index = 0;
        for(int j = 0; j < indexList.size(); j+=3){
            Vec3 A = vertList[indexList[j]];
            Vec3 B = vertList[indexList[j+1]];
            Vec3 C = vertList[indexList[j+2]];

            //get midpoints of AB, AC, BC, move them out along the radius
            Vec3 AB = getMidpoint(A,B);
            AB = origin-AB;
            AB = AB.normalized();
            AB = radius*AB;
            AB = origin - AB;

            Vec3 AC = getMidpoint(A,C);
            AC = origin-AC;
            AC = AC.normalized();
            AC = radius*AC;
            AC = origin - AC;

            Vec3 BC = getMidpoint(B,C);
            BC = origin-BC;
            BC = BC.normalized();
            BC = radius*BC;
            BC = origin - BC;

            tempVertList.push_back(A);  //0
            tempVertList.push_back(B);  //1
            tempVertList.push_back(C);  //2

            tempVertList.push_back(AB); //3
            tempVertList.push_back(AC); //4
            tempVertList.push_back(BC); //5

            //add (A,AB,AC),(B,BC,AB),(C,AC,BC),(AB,BC,AC) to temp indexList
            tempIndexList.push_back(index);
            tempIndexList.push_back(index+3);
            tempIndexList.push_back(index+4);

            tempIndexList.push_back(index+1);
            tempIndexList.push_back(index+3);
            tempIndexList.push_back(index+5);

            tempIndexList.push_back(index+2);
            tempIndexList.push_back(index+4);
            tempIndexList.push_back(index+5);

            tempIndexList.push_back(index+3);
            tempIndexList.push_back(index+4);
            tempIndexList.push_back(index+5);

            index+=6;

        }

        //replace lists with temps
        vertList = tempVertList;
        indexList = tempIndexList;
    }

    renderMesh.loadVertices(vertList, indexList);

    getNormals(vertList,indexList,normList);

    renderMesh.loadNormals(normList);

    /// Load textures (assumes texcoords)
    renderMesh.loadTextures(texFile);

    /// Load texture coordinates (assumes textures)
    /// TODO: add texcoords
    for(int i = 0; i < vertList.size(); i++){
        Vec2 uv = sphereMap(vertList[i], origin, radius);
        tCoordList.push_back(uv);
    }
    renderMesh.loadTexCoords(tCoordList);

}

/// Loads a cylinder into renderMesh via vertList, indexList, and normList, with the given origin point, radius, height, & subvidision level
void loadCylinder(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<Vec2> &tCoordList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, const char * texFile, Vec3 origin, float radius, float height, int sublevel){
    if(sublevel < 4){
        sublevel = 4;
    }
    int numVerts = (2*sublevel)+2;
    float vertices[numVerts][3];
    float angle = (2*pi)/sublevel;
    float the = 0.0f;

    //add top & bottom verts
    vertices[0][0] = 0.0f;
    vertices[0][1] = 0.0f;
    vertices[0][2] = 0.0f;
    vertices[1][0] = 0.0f;
    vertices[1][1] = height;
    vertices[1][2] = 0.0f;

    for(int i = 2; i < numVerts; i+=2){
        //std::cout << i << std::endl;
        float x = radius*cos(the);
        float z = radius*sin(the);
        vertices[i][0] = x;
        vertices[i][1] = 0.0f;
        vertices[i][2] = z;
        vertices[i+1][0] = x;
        vertices[i+1][1] = height;
        vertices[i+1][2] = z;
        the+=angle;
    }

    for(int i = 0; i < numVerts; i++){
        vertList.push_back(Vec3(vertices[i][0],vertices[i][1],vertices[i][2]));
    }

    //add top triangle indices
    for(int i = 4; i < numVerts-1; i+=2){
        indexList.push_back(0);
        indexList.push_back(i-2);
        indexList.push_back(i);
    }
    indexList.push_back(0);
    indexList.push_back(2);
    indexList.push_back(2*sublevel);

    //add bottom triangle indices
    for(int i = 5; i < numVerts; i+=2){
        indexList.push_back(1);
        indexList.push_back(i-2);
        indexList.push_back(i);
    }
    indexList.push_back(1);
    indexList.push_back(3);
    indexList.push_back((2*sublevel)+1);

    for(int i = 2; i < numVerts-2; i++){
        indexList.push_back(i);
        indexList.push_back(i+1);
        indexList.push_back(i+2);
    }
    indexList.push_back(numVerts-2);
    indexList.push_back(numVerts-1);
    indexList.push_back(2);

    indexList.push_back(2);
    indexList.push_back(3);
    indexList.push_back(numVerts-1);

    renderMesh.loadVertices(vertList, indexList);

    getNormals(vertList,indexList,normList);

    renderMesh.loadNormals(normList);

    /// Load textures (assumes texcoords)
    renderMesh.loadTextures(texFile);

    /// Load texture coordinates (assumes textures)
    /// TODO: add texcoords
    renderMesh.loadTexCoords(tCoordList);
}

/// Loads a .obj file from the given filepath, transfers position, texture, index, & normal data into vectors, loads vectors into renderMesh
/// Based off code from: http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/
bool loadObj(const char * path, Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<Vec2> &tCoordList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, const char * texFile){

    FILE * file = fopen(path,"r");

    if(file == NULL){
        printf("File not found.");
        return false;
    }

    while(1){
        char lineHeader[128];

        int res = fscanf(file, "%s", lineHeader);
        if(res == EOF){
            break;
        }
        //issue: generates bizarre wireframe
        ///TODO: add in vn, vt, and bools to track them

        if(strcmp(lineHeader, "v") == 0){
            Vec3 vertex;
            fscanf(file, "%f %f %f\n",&vertex[0],&vertex[1],&vertex[2]);
            vertList.push_back(vertex);
        }else if(strcmp(lineHeader, "f") == 0){
            unsigned int v1,v2,v3;
            fscanf(file, "%d %d %d\n",&v1,&v2,&v3);
            indexList.push_back(v1);
            indexList.push_back(v2);
            indexList.push_back(v3);
        }
    }

    renderMesh.loadVertices(vertList, indexList);
    renderMesh.loadNormals(normList);

    /// Load textures (assumes texcoords)
    renderMesh.loadTextures(texFile);

    /// Load texture coordinates (assumes textures)
    renderMesh.loadTexCoords(tCoordList);

    return true;
}

/// Writes a .obj file from the current renderMesh to the build directory
void writeObj(std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList){

    std::ofstream ofs ("output.obj",std::ofstream::out);

    for(int i = 0; i < vertList.size(); i++){
        ofs << "v " << vertList[i][0] << " " << vertList[i][1] << " " << vertList[i][2] << std::endl;
    }

    ofs << std::endl;

    for(int i = 0; i < indexList.size(); i+= 3){
        ofs << "f " << indexList[i] << " " << indexList[i+1] << " " << indexList[i+2] << std::endl;
    }

    ofs.close();

}

int main() {

    Application app;
    ImguiRenderer imrenderer;
    Mesh renderMesh;

    /// Example rendering a mesh
    /// Call to compile shaders
    renderMesh.init();

    /// Initialize vector/textures/index/normal lists
    std::vector<Vec3> vertList;
    std::vector<Vec2> tCoordList;
    std::vector<unsigned int> indexList;
    std::vector<Vec3> normList;

    Vec3 origin = Vec3(0.0f,0.0f,0.0f);
    float size = 1.0f;
    float height = 0.1f;
    const char * objFile = "bunny.obj";
    const char * texFile = "earth.png";

    /// Load cube vertices, indices, and normals (pending)
    //loadCube(renderMesh, vertList, tCoordList, indexList, normList, texFile, origin, size);

    /// Load sphere vertices, indices, and normals (pending)
    loadIcoSphere(renderMesh, vertList, tCoordList, indexList, normList, texFile, origin, size, 4);

    /// Load cylinder vertices, indices, and normals (pending)
    //loadCylinder(renderMesh, vertList, tCoordList, indexList, normList, texFile, origin, size, height, 100);

    /// Load mesh from .obj filepath
    //loadObj(objFile, renderMesh, vertList, tCoordList, indexList, normList);

    writeObj(vertList,indexList);

    /// Create main window, set callback function
    auto &window1 = app.create_window([&](Window &window){
        int width, height;
        std::tie(width, height) = window.get_size();

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glClearColor(0.0f, 0.0f, 0.0f, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        /// Wireframe rendering, might be helpful when debugging your mesh generation
        //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

        float ratio = width / (float) height;
        Mat4x4 modelTransform = Mat4x4::Identity();
        Mat4x4 model = modelTransform.matrix();
        Mat4x4 projection = OpenGP::perspective(70.0f, ratio, 0.1f, 10.0f);

        //camera movement
        float time = .5f * (float)glfwGetTime();
        Vec3 cam_pos(2*cos(time), 2.0, 2*sin(time));
        Vec3 cam_look(0.0f, 0.0f, 0.0f);
        Vec3 cam_up(0.0f, 1.0f, 0.0f);
        Mat4x4 view = OpenGP::lookAt(cam_pos, cam_look, cam_up);

        renderMesh.draw(model, view, projection);
    });
    window1.set_title("Assignment 2");

    /*
    /// Create window for IMGUI, set callback function
    auto &window2 = app.create_window([&](Window &window){
        int width, height;
        std::tie(width, height) = window.get_size();
        imrenderer.begin_frame(width, height);
        ImGui::BeginMainMenuBar();
        ImGui::MenuItem("File");
        ImGui::MenuItem("Edit");
        ImGui::MenuItem("View");
        ImGui::MenuItem("Help");
        ImGui::EndMainMenuBar();
        ImGui::Begin("Test Window 1");
        ImGui::Text("This is a test imgui window");
        ImGui::End();
        glClearColor(0.15f, 0.15f, 0.15f, 1);
        glClear(GL_COLOR_BUFFER_BIT);
        imrenderer.end_frame();
    });
    window2.set_title("imgui Test");
    */

    return app.run();
}
