# ============================================================================
#  simlab — jour 7
# ============================================================================
#  Cibles :
#    simcore            bibliothèque numérique, AUCUNE dépendance OpenGL
#    simlab_tests       tests doctest (doctest.h est fourni dans tests/)
#    advection_demo, projection_demo, smoke_demo, pic_transfer_demo
#                       démos console -> CSV
#    simlab             fumée interactive (jour 5)
#    simlab_particles   particules PIC interactives (jour 7)
#    (les deux applications : option SIMLAB_BUILD_GUI)
#
#  GLFW et ImGui sont téléchargés et compilés par CMake (FetchContent) à la
#  première configuration. glad (OpenGL 4.3 Core) est fourni dans external/glad.
# ============================================================================
cmake_minimum_required(VERSION 3.16)
project(simlab C CXX)          # C : pour glad.c

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- Cœur numérique : aucune dépendance OpenGL ---------------------------------
add_library(simcore
  src/grid/Field2D.cpp
  src/grid/MacGrid2D.cpp
  src/numerics/Advection.cpp
  src/numerics/MacOperators.cpp
  src/numerics/ConjugateGradient.cpp
  src/numerics/Projection.cpp
  src/numerics/Diffusion.cpp               # J5
  src/numerics/Forces.cpp                  # J5
  src/diagnostics/FluidDiagnostics.cpp     # J5
  src/solvers/StableFluidsSolver.cpp       # J5
  src/scenes/Seeding.cpp                   # J7
  src/numerics/MacTransfer.cpp             # J7
  src/diagnostics/ParticleDiagnostics.cpp  # J7
  src/solvers/FlipSolver.cpp)              # J7
target_include_directories(simcore PUBLIC src)
if(MSVC)
  target_compile_options(simcore PUBLIC /W4 /utf-8)
else()
  target_compile_options(simcore PUBLIC -Wall -Wextra)
endif()

# --- Démos console ----------------------------------------------------------------
foreach(demo advection_demo projection_demo smoke_demo pic_transfer_demo)   # pic_transfer_demo : J7
  add_executable(${demo} src/app/${demo}.cpp)
  target_link_libraries(${demo} PRIVATE simcore)
endforeach()

# --- Tests unitaires ------------------------------------------------------------
enable_testing()
add_executable(simlab_tests
  tests/test_main.cpp
  tests/test_grid.cpp
  tests/test_advection.cpp
  tests/test_projection.cpp
  tests/test_smoke.cpp                     # J5
  tests/test_pic.cpp)                      # J7
target_include_directories(simlab_tests PRIVATE tests)
target_link_libraries(simlab_tests PRIVATE simcore)
add_test(NAME simlab_tests COMMAND simlab_tests)

# --- Applications interactives ----------------------------------------------------
option(SIMLAB_BUILD_GUI "Applications interactives (GLFW + glad + ImGui)" ON)
if(SIMLAB_BUILD_GUI)
  include(FetchContent)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)    # pour les CMake récents
  set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_WAYLAND  OFF CACHE BOOL "" FORCE)   # Linux seulement : X11 suffit
  FetchContent_Declare(glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
    GIT_SHALLOW TRUE)
  FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.5
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(glfw imgui)   # ImGui n'a pas de CMakeLists : simple téléchargement

  find_package(OpenGL REQUIRED)

  add_library(glad STATIC external/glad/src/glad.c)
  target_include_directories(glad PUBLIC external/glad/include)

  add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
  target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
  target_link_libraries(imgui PUBLIC glfw)

  # Tout le code OpenGL du projet vit ici (src/render) et dans les main.
  add_library(render_gl STATIC
    src/render/FieldRenderer.cpp
    src/render/ParticleRenderer.cpp)       # J7
  target_link_libraries(render_gl PUBLIC simcore glad OpenGL::GL)

  foreach(app simlab simlab_particles)
    if(app STREQUAL "simlab")
      add_executable(${app} src/app/main.cpp)
    else()
      add_executable(${app} src/app/particles_main.cpp)   # J7
    endif()
    target_link_libraries(${app} PRIVATE render_gl imgui glfw)
    # Chemin absolu des shaders : l'exécutable les trouve quel que soit le dossier courant.
    target_compile_definitions(${app} PRIVATE SIMLAB_SHADER_DIR="${CMAKE_CURRENT_SOURCE_DIR}/shaders/")
  endforeach()
endif()
