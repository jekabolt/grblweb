var svgdiv = Snap("#svgdiv");
// var contentGroup = svgdiv.group();
var contentGroup = svgdiv.group();

var toolPathsGroup = svgdiv.g();
var selectionGroup = svgdiv.g();
var tabsGroup = svgdiv.g();
var combinedGeometryGroup = svgdiv.g();

var miscViewModel = new MiscViewModel();
var svgViewModel = new SvgViewModel();
var toolModel = new ToolModel();

var materialViewModel = new MaterialViewModel();
var selectionViewModel = new SelectionViewModel(svgViewModel, materialViewModel, selectionGroup);

var gcodeConversionViewModel = new GcodeConversionViewModel(options, miscViewModel, materialViewModel, toolModel,
    operationsViewModel, tabsViewModel);

var tabsViewModel = new TabsViewModel(
    miscViewModel, options, svgViewModel, materialViewModel, selectionViewModel, tabsGroup,
    function() { gcodeConversionViewModel.generateGcode(); });

var operationsViewModel = new OperationsViewModel(
    miscViewModel, options, svgViewModel, materialViewModel, selectionViewModel, toolModel, combinedGeometryGroup, toolPathsGroup,
    function() { gcodeConversionViewModel.generateGcode(); });



window.onload = function() {
    var imageSelect = document.getElementById("imageSelect"),
        imageInput = document.getElementById("imageInput");

    imageSelect.addEventListener("click", function(e) {
        imageInput.click();
        e.preventDefault();
    }, false);

    imageInput.addEventListener("change", function(e) {
        handleFiles(this.files);
    }, false);

    document.getElementById("save").addEventListener("click", function(e) {
        e.target.download = "potrace" + (new Date()).toLocaleTimeString() + ".svg";
        e.target.href = "data:image/svg+xml;," + Potrace.getSVG(1);
    }, false);

    document.getElementById("svgdiv").onclick = function(e) {
        var element = Snap.getElementByPoint(457, 515);
        console.log(e.pageX, e.pageY);
        console.log(element);

        // if (element != null) {
        //     operationsViewModel.clickOnSvg(element) || tabsViewModel.clickOnSvg(element) || selectionViewModel.clickOnSvg(element);
        //     if (selectionViewModel.selNumSelected() > 0) {
        //         tutorial(3, 'Click "Create Operation" after you have finished selecting objects.');
        //     }
        // }
    }

    document.getElementById("act").onclick = function(e) {
        // console.log(contentGroup.node.childNodes);
        trap.kek = "ttrrraaaappp"
        console.log(trap);
    }



};

function handleFiles(files) {
    Potrace.loadImageFromFile(files[0]);
    Potrace.process(function() {
        var reader = new FileReader();
        console.log(reader);
        reader.onload = function(e) {
            loadSvg(Potrace.getSVG(1));
        };

        var blob = new Blob([Potrace.getSVG(1)]);
        reader.readAsText(blob);
        displaySVG(1);
    });
}

function loadSvg(content) {
    svg = Snap.parse(content);
    contentGroup.append(svg);
    bbox = svgdiv.getBBox();
    document.getElementById("svgdiv").setAttribute("viewBox", (bbox.x - 2) + " " + (bbox.y - 2) + " " + (bbox.w + 4) + " " + (bbox.h + 4));
    // updateSvgSize();
}

function updateSvgSize() {
    bbox = mainSvg.getBBox();
    $("#MainSvg").attr({
        width: $("#MainSvgDiv").width(),
        height: Math.max(10, $(window).height() - 120),
        preserveAspectRatio: 'xMinYMin meet',
    });
    // attr() messes viewBox up
    $("#MainSvg").get(0).setAttribute("viewBox", (bbox.x - 2) + " " + (bbox.y - 2) + " " + (bbox.w + 4) + " " + (bbox.h + 4));
}

function displaySVG(size, type) {
    var svgdiv = document.getElementById('svgdiv');
    svgdiv.style.display = 'inline-block';
    svgdiv.innerHTML = "<p>Result:</p>" + Potrace.getSVG(size, type);
    save(Potrace.getSVG(size, type), "vector", ".svg")
}

function save(svg, name, type) {
    var dlbtn = document.getElementById("dlbtn");
    var file = new Blob([svg], { type: type });
    dlbtn.href = URL.createObjectURL(file);
    dlbtn.download = name;
}

function MiscViewModel() {
    var self = this;
    // self.debug = ko.observable(options.debug);

    // self.saveGistDescription = ko.observable("jscut settings");
    // self.savedGistUrl = ko.observable("");
    // self.savedGistLaunchUrl = ko.observable("");
    // self.localStorageSettings = ko.observableArray([]);
    // self.loadedCamCpp = ko.observable(false);
    // self.camCppError = ko.observable("");
}

var options = {
    camCppPaths: ['js', 'http://api.jscut.org/js'],
    enableGoogleDrive: true,
    enableDropbox: true,
    profile: false,
    debug: false,
    preloadInBrowser: 'preload.jscut',
};