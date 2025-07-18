let selectedBlockElement = null;
let floatingInfoBox = null;


function initAllChromosomesVisualization() {
    var chartContainer = document.getElementById("chart");

    // Очищаем контейнер перед отрисовкой!
    chartContainer.innerHTML = "";

    const drawnChromosomes = new Set();

    for (var id in references_by_id) {
        var chrName = references_by_id[id];

        if (drawnChromosomes.has(chrName)) {
            console.warn(`Хромосома ${chrName} уже была отрисована! Пропускаем повтор.`);
            continue;
        }

        drawnChromosomes.add(chrName);

        var chrDiv = document.createElement("div");
        chrDiv.className = "chromosome_chart";
        chrDiv.id = "chart_" + chrName;
        chartContainer.appendChild(chrDiv);

        drawChromosome(chrName, chrDiv.id);
    }
}



function drawChromosome(chrName, containerId) {
    var contigs = contig_data[chrName];
    var chrLength = chromosomes_len[chrName];
    var svgWidth = 1000;
    var laneHeight = 22;
    var marginTop = 30;
    var marginLeft = 60;
    var assemblyLabels = Object.keys(contigs);
    var svgHeight = laneHeight * assemblyLabels.length + marginTop + 20;
    var colorScale = d3.scale.category10();

    var svg = d3.select("#" + containerId)
        .append("svg")
        .attr("width", svgWidth)
        .attr("height", svgHeight);

    // pattern for centromere
    svg.append("defs").append("pattern")
        .attr("id", "centromerePattern")
        .attr("patternUnits", "userSpaceOnUse")
        .attr("width", 6)
        .attr("height", 6)
        .append("path")
        .attr("d", "M0,0 l6,6")
        .attr("stroke", "#888")
        .attr("stroke-width", 1);

    // X-axis scale
    var xScale = d3.scale.linear()
        .domain([0, chrLength])
        .range([marginLeft, svgWidth - 40]);

    // X-axis
    var xAxis = d3.svg.axis()
        .scale(xScale)
        .orient('bottom')
        .ticks(6)
        .tickFormat(formatLongNumbers);

    var xAxisGroup = svg.append('g')
        .attr('class', 'x axis')
        .attr('transform', 'translate(0,' + (svgHeight - 20) + ')')
        .call(xAxis);

    var ticks = xAxisGroup.selectAll(".tick text");

    if (!ticks.empty()) {
        var lastTick = ticks[0][ticks[0].length - 1];
        d3.select(lastTick).text(formatLongNumbers(chrLength));
    }

    xAxisGroup.selectAll("path, line")
        .style("fill", "none")
        .style("stroke", "#000")
        .style("stroke-width", "0.5px");

    assemblyLabels.forEach(function (label, idx) {
        var contigList = contigs[label];
        var y0 = marginTop + laneHeight * idx;

        // Имя дорожки (assembly)
        svg.append("text")
            .attr("x", 0)
            .attr("y", marginTop + laneHeight * idx + 15)
            .text(label)
            .attr("font-size", "12px")
            .attr("fill", colorScale(idx));

        // Прямоугольники (блоки)
        svg.selectAll("rect_" + label)
            .data(contigList)
            .enter()
            .append("rect")
            .attr("x", function (d) { return xScale(d.start); })
            .attr("y", marginTop + laneHeight * idx)
            .attr("width", function (d) { return Math.max(1, xScale(d.end) - xScale(d.start)); })
            .attr("height", laneHeight - 4)
            .attr("fill", function (d) {
                if (d.name && d.name.toLowerCase().includes("centromere")) {
                    return "url(#centromerePattern)";
                }
                return colorScale(idx);
            })
            .attr("class", function (d) {
                let laneClass = "lane_" + label;
                let centromereClass = "";
                if (d.name && d.name.toLowerCase().includes("centromere")) {
                    centromereClass = " centromere_block";
                }
                if (d.misassembledEnds)
                    return "block end " + (d.objClass || "") + " " + laneClass + centromereClass;
                if (!d.marks || d.contig_type)
                    return "block mainItem " + (d.objClass || "") + " " + laneClass + centromereClass;
                return "block " + laneClass + centromereClass;
            })
            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)
            .attr("data-name", d => d.name)
            .on("click", function (d, i) {
                console.log("[Клик по rect]:", d);
                handleBlockClick(d, this);
            })
            .append("title")
            .text(function (d) {
                return d.name + ": " + d.start + " - " + d.end;
            });

        // ==== Генерируем triangles массив ====
        const triangles = [];
        contigList.forEach(d => {
            // Только если mis_ends не пустое
            if (d.mis_ends && d.mis_ends.trim() !== "") {
                d.mis_ends.split(';').forEach(end => {
                    end = end.trim();
                    if (end === "L" || end === "R") {
                        triangles.push({
                            ...d,
                            mis_end_side: end
                        });
                        console.log(`[${chrName} | ${label}] ТРЕУГОЛЬНИК создан для`, d.name, 'side:', end, d);
                    }
                });
            }
        });


        svg.selectAll(".triangle_" + label)
            .data(triangles)
            .enter()
            .append("path")
            .attr("class", "block end misassembled mis_triangle triangle_" + label)
            .attr("d", function (d) {
                const size = 7;
                const xStart = xScale(d.start);
                const xEnd = xScale(d.end);
                const yTop = marginTop + laneHeight * idx;
                const yMid = yTop + (laneHeight - 4) / 2;

                if (d.mis_end_side === "L") {
                    return `M${xStart},${yMid - size} L${xStart},${yMid + size} L${xStart + size},${yMid} Z`;
                } else if (d.mis_end_side === "R") {
                    return `M${xEnd},${yMid - size} L${xEnd},${yMid + size} L${xEnd - size},${yMid} Z`;
                }
                return null;
            })
            .attr("fill", "#cc2222")
            .attr("stroke", "#222")
            .attr("stroke-width", 1)
    });  // ВОТ ТУТ закрой forEach!

    // ---- Текст названия хромосомы ----
    svg.append("text")
        .attr("x", 5)
        .attr("y", 18)
        .text(chrName)
        .attr("font-size", "16px")
        .attr("fill", "#333");
}


let selectedBlockName = null;

function handleBlockClick(blockData, element) {
    const name = blockData.name;
    const clicked = d3.select(element);

    d3.selectAll('rect[data-name]')
        .attr('stroke-width', 0.5)
        .attr('stroke', '#000')
        .attr('stroke-dasharray', null);

    if (selectedBlockName === name) {
        selectedBlockName = null;
        clearInfoPanel();
        if (floatingInfoBox) { floatingInfoBox.remove(); floatingInfoBox = null; }
        return;
    }


    d3.selectAll(`rect[data-name="${name}"]`)
        .filter(function () { return this !== element; })
        .attr('stroke', '#f00')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '8,4');


    clicked
        .attr('stroke-width', 2.5)
        .attr('stroke', '#f00')
        .attr('stroke-dasharray', null);

    selectedBlockName = name;


    updateInfoPanel(blockData);

    if (floatingInfoBox) {
        floatingInfoBox.remove();
        floatingInfoBox = null;
    }

    const svg = d3.select(element.ownerSVGElement);
    const x = +clicked.attr("x");
    const y = +clicked.attr("y");
    const lines = [
        `${blockData.name || "Contig"}, Start: ${formatLongNumbers(blockData.start)}, End: ${formatLongNumbers(blockData.end)}`
    ];

    const tempSvg = d3.select("body").append("svg").attr("visibility", "hidden");
    let maxWidth = 0;
    lines.forEach(line => {
        const tempText = tempSvg.append("text").attr("font-size", "12px").text(line);
        const length = tempText.node().getComputedTextLength();
        if (length > maxWidth) maxWidth = length;
        tempText.remove();
    });
    tempSvg.remove();
    const boxWidth = maxWidth + 20;
    const lineHeight = 14;
    const boxHeight = lines.length * lineHeight + 8;

    const svgWidth = +svg.node().getAttribute("width");
    let adjustedX = x;
    if (x + boxWidth > svgWidth - 10) {
        adjustedX = svgWidth - boxWidth - 10;
    }

    floatingInfoBox = svg.append("g").attr("class", "floating-info");
    floatingInfoBox.append("rect")
        .attr("x", adjustedX)
        .attr("y", y - boxHeight - 5)
        .attr("width", boxWidth)
        .attr("height", boxHeight)
        .attr("fill", "#fff")
        .attr("stroke", "#333")
        .attr("stroke-width", 1)
        .attr("rx", 4).attr("ry", 4);

    lines.forEach((line, i) => {
        floatingInfoBox.append("text")
            .attr("x", adjustedX + 6)
            .attr("y", y - boxHeight + lineHeight * i + 8)
            .attr("font-size", "12px")
            .attr("fill", "#000")
            .text(line);
    });
}



function formatLongNumbers(d) {
    if (d >= 1000000) {
        return (d / 1000000).toFixed(2).replace(/\.0+$/, '') + " Mbp";
    } else if (d >= 1000) {
        return (d / 1000).toFixed(1).replace(/\.0$/, '') + " Kbp";
    } else {
        return d + " bp";
    }
}

function updateInfoPanel(blockData) {
    var infoPanel = document.getElementById("block_info");
    if (infoPanel) {
        let chrSize = chromosomes_len[blockData.name];
        let sizeStr = chrSize ? " (" + formatLongNumbers(chrSize) + ")" : "";

        infoPanel.innerHTML =
            `<b>${blockData.name || "?"}${sizeStr}</b><br/>
            Start: ${blockData.start || "?"}<br/>
            End: ${blockData.end || "?"}<br/>
            Contig Type: ${blockData.contig_type || "Unknown"}<br/>
            Is Best: ${blockData.is_best || "False"}<br/>
            Misassemblies: ${blockData.misassemblies || "None"}`;
    } else {
        console.warn("info_panel not found in DOM");
    }
}

function clearInfoPanel() {
    var infoPanel = document.getElementById("block_info");
    if (infoPanel) {
        infoPanel.innerHTML = "";
    }
}


document.addEventListener("DOMContentLoaded", function () {
    initAllChromosomesVisualization();
    window.items = [];
    for (var chr in contig_data) {
        var contigSets = contig_data[chr];
        for (var label in contigSets) {
            var blocks = contigSets[label];
            window.items.push(...blocks);

        }
    }

    var itemsLayer = d3.select('body')
        .append('svg')
        .attr('id', 'items');
    window.itemsContainer = itemsLayer.append('g');

    if (typeof setupInterface !== 'undefined') {
        setupInterface();
    } else {
        console.warn("setupInterface() не найден");
    }

    console.log("Общее количество items:", items.length);
});