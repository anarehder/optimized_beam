render(convexity = 10) difference() {
import("DISP_TESTE.stl", convexity=10);
translate([-30, -15, -1.5]) {
    import(  "/home/ana/Algoritimo_Genetico_Opt_Beam/build/DISP_IDEAL.stl");
    };
    } 