const pptxgen = require('pptxgenjs');
const html2pptx = require('./scripts/html2pptx');

(async () => {
    try {
        console.log("=== ANTHROPIC SKILLS DEMO ===");
        console.log("Initializing local browser connection...");
        
        const pres = new pptxgen();
        pres.layout = 'LAYOUT_16x9';

        console.log("Converting 'slide1.html' (Cover) to PowerPoint slide...");
        await html2pptx('slide1.html', pres);

        console.log("Converting 'slide2.html' (Content) to PowerPoint slide...");
        await html2pptx('slide2.html', pres);
        
        console.log("Saving 'My_Self_Intro.pptx'...");
        await pres.writeFile('My_Self_Intro.pptx');
        console.log("✅ SUCCESS: Presentation generated successfully!");
    } catch (err) {
        console.error("❌ ERROR:", err);
    }
})();