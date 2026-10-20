#@ String inDir
#@ String outFile
// Fiji script that produced coloc2_golden.csv: runs Coloc 2 (Colocalisation_Analysis)
// on every <case>_ch1.tif / <case>_ch2.tif pair in inDir, once per threshold-regression
// implementation, and writes every reported value as  case,implementation,name,value.
// Settings: no ROI/mask, PSF 3, 100 Costes randomisations, Manders + Costes test on,
// Li / Spearman / Kendall / scatterplot off.
// Run headless (see makeColocGolden.m):
//   fiji-windows-x64.exe --headless --run coloc2Golden.groovy "inDir='...',outFile='...'"
import ij.IJ
import sc.fiji.coloc.Coloc_2
import sc.fiji.coloc.algorithms.AutoThresholdRegression

def out = new File(outFile)
def lines = ["case,implementation,name,value"]
def meta = [
    "ij=" + IJ.getFullVersion(),
    "coloc2jar=" + Coloc_2.class.getProtectionDomain().getCodeSource().getLocation(),
    "java=" + System.getProperty("java.version")]
try {
    def files = new File(inDir).listFiles().findAll { it.name.endsWith("_ch1.tif") }.sort { it.name }
    for (f1 in files) {
        def caseName = f1.name.replace("_ch1.tif", "")
        def f2 = new File(inDir, caseName + "_ch2.tif")
        AutoThresholdRegression.Implementation.values().eachWithIndex { impl, iImpl ->
            def imp1 = IJ.openImage(f1.path)
            def imp2 = IJ.openImage(f2.path)
            def c = new Coloc_2()
            // (imp1, imp2, mask index 0 = <None>, regression index, autoSavePdf, displayImages,
            //  displayShuffledCostes, Li ch1, Li ch2, Li ICQ, Spearman, Manders, Kendall,
            //  scatterplot, Costes test, PSF, Costes randomisations)
            def ok = c.initializeSettings(imp1, imp2, 0, iImpl, false, false, false, false, false,
                false, false, true, false, false, true, 3, 100)
            if (!ok) {
                lines << [caseName, impl.toString(), "error", "initializeSettings failed"].join(",")
                return
            }
            def res = c.colocalise(c.img1, c.img2, null, null, null)
            for (v in res.values()) {
                def val = v.isNumber ? String.valueOf(v.number) : ("\"" + v.value + "\"")
                lines << [caseName, impl.toString(), "\"" + v.name + "\"", val].join(",")
            }
            for (w in res.warnings()) {
                lines << [caseName, impl.toString(), "\"warning\"", "\"" + w.getShortMessage() + "\""].join(",")
            }
        }
    }
} catch (Throwable e) {
    def sw = new StringWriter()
    e.printStackTrace(new PrintWriter(sw))
    lines << "ERROR," + "\"" + sw.toString().replace("\"", "'").replace("\n", " | ") + "\""
}
out.text = "# " + meta.join("; ") + "\n" + lines.join("\n") + "\n"
