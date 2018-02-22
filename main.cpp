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

/// Loads a cube into renderMesh via vertList, indexList, and normList, with the given origin point & side length
void loadCube(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, Vec3 origin, float side){
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
}

/// Loads a UVsphere into renderMesh via vertList, indexList, and normList, with the given origin point, radius, & level of subdivision
void loadUVSphere(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, Vec3 origin, float radius, int sublevel){
   // float vertices[][];
    //from the origin, make a vector with radius length & (0,0,0) direction
    //do subdivision times:
        //do subdivision/2 times:
}

/// Loads an icosphere into renderMesh via vertList, indexList, and normList, with the given origin point, radius, & level of subdivision (pending)
/// based off pseudocode from: https://www.csee.umbc.edu/~squire/reference/polyhedra.shtml#icosahedron
void loadIcoSphere(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, Vec3 origin, float radius, int sublevel){
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

    renderMesh.loadVertices(vertList, indexList);

    getNormals(vertList,indexList,normList);

    renderMesh.loadNormals(normList);

}

/// Loads a cylinder into renderMesh via vertList, indexList, and normList, with the given origin point, radius, height, & subvidision level
void loadCylinder(Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList, Vec3 origin, float radius, float height, int sublevel){
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
}

/// Loads a .obj file from the given filepath, transfers position, index, & normal data into vectors, loads vectors into renderMesh
/// Based off code from: https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Load_OBJ
bool loadObjIfstream(const char * path, Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<unsigned int> &indexList){
    std::ifstream in(path, std::ios::in);
    if(!in){
        printf("File not found.");
        return false;
    }

    std::string line;
    while(getline(in,line)){
        if(line.substr(0,2) == "v"){
            std::istringstream s(line.substr(2));
            Vec3 v; s >> v[0]; s >> v[1]; s >> v[2];
            vertList.push_back(v);
        }else if(line.substr(0,2) == "f"){
            std::istringstream s(line.substr(2));
            unsigned int a,b,c;
            s >> a; s >> b; s >> c;
            indexList.push_back(a);
            indexList.push_back(b);
            indexList.push_back(c);
        }
    }

    renderMesh.loadVertices(vertList, indexList);

}

/// Loads a .obj file from the given filepath, transfers position, texture, index, & normal data into vectors, loads vectors into renderMesh
/// Based off code from: http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/
bool loadObjFopen(const char * path, Mesh &renderMesh, std::vector<Vec3> &vertList, std::vector<Vec2> &tCoordList, std::vector<unsigned int> &indexList, std::vector<Vec3> &normList){

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

        if(strcmp(lineHeader, "v") == 0){
            Vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex[0],&vertex[1],&vertex[2]);
            vertList.push_back(vertex);
        }/*else if(strcmp(lineHeader, "vt") == 0){
            Vec2 uv;
            fscanf(file, "%f %f\n", &uv[0], &uv[1]);
            tCoordList.push_back(uv);
        }else if(strcmp(lineHeader, "vn") == 0){
            Vec3 normal;
            fscanf(file, "%f %f %f\n", &normal[0], &normal[1], &normal[2]);
            normList.push_back(normal);
        }*/else if(strcmp(lineHeader, "f") == 0){
            unsigned int vertexIndex[3];
            fscanf(file, "%d/%d/%d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);  //update this to work with %d/%d/%d
            indexList.push_back(vertexIndex[0]);
            indexList.push_back(vertexIndex[1]);
            indexList.push_back(vertexIndex[2]);
        }
    }

    renderMesh.loadVertices(vertList, indexList);
    renderMesh.loadNormals(normList);
}

int main() {

    Application app;
    ImguiRenderer imrenderer;
    Mesh renderMesh;

    /// Example rendering a mesh
    /// Call to compile shaders
    renderMesh.init();

    /// TODO: write functions to create different objects:
        /// -sphere
        /// -cylinder
        /// -read in .obj file

    /// Initialize vector/textures/index/normal lists
    std::vector<Vec3> vertList;
    std::vector<Vec2> tCoordList;
    std::vector<unsigned int> indexList;
    std::vector<Vec3> normList;

    Vec3 origin = Vec3(0.0f,0.0f,0.0f);
    float size = 1.0f;
    float height = 2.0f;

    /// Load cube vertices, indices, and normals (pending)
    //loadCube(renderMesh, vertList, indexList, normList, origin, size);

    /// Load sphere vertices, indices, and normals (pending)
    //loadIcoSphere(renderMesh, vertList, indexList, normList, origin, size, 1);

    /// Load cylinder vertices, indices, and normals (pending)
    //loadCylinder(renderMesh, vertList, indexList, normList, origin, size, height, 12);

    /// Load mesh from .obj filepath
    const char * path = "bunny.obj";
    //loadObjFopen(path, renderMesh, vertList, tCoordList, indexList, normList);
    //loadObjIfstream(path, renderMesh, vertList, indexList);

    /// TODO: get textures working

    /// Load textures (assumes texcoords)
    renderMesh.loadTextures("earth.png");

    /// Load texture coordinates (assumes textures)
    tCoordList.push_back(Vec2(0,0));
    tCoordList.push_back(Vec2(1,0));
    tCoordList.push_back(Vec2(1,1));
    tCoordList.push_back(Vec2(0,1));
    renderMesh.loadTexCoords(tCoordList);

    /// Create main window, set callback function
    auto &window1 = app.create_window([&](Window &window){
        int width, height;
        std::tie(width, height) = window.get_size();

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glClearColor(0.0f, 0.0f, 0.0f, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        /// Wireframe rendering, might be helpful when debugging your mesh generation
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

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
