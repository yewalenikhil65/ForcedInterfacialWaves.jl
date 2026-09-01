/* Pair the existing Julia and MATLAB code blocks on theory pages.
 * Executable Julia examples may generate a <pre> output block or an image
 * before the MATLAB block; those nodes stay with the Julia column. */
(function () {
    function isLanguageBlock(node, language) {
        return node && node.matches("pre") &&
            node.querySelector("code.language-" + language);
    }

    function isGeneratedJuliaOutput(node) {
        return node && (
            node.matches("pre.documenter-example-output") ||
            node.matches("img") ||
            node.matches("svg") ||
            node.classList.contains("documenter-example-output")
        );
    }

    function isBoundary(node) {
        return node && node.matches("h1, h2, h3, h4, h5, h6, hr, p, ul, ol, table, figure");
    }

    function pairExamples() {
        document.querySelectorAll("pre").forEach(function (juliaBlock) {
            if (!isLanguageBlock(juliaBlock, "julia") ||
                juliaBlock.closest(".code-comparison")) {
                return;
            }

            var matlabBlock = null;
            var generatedOutput = [];
            var node = juliaBlock.nextElementSibling;

            while (node) {
                if (isLanguageBlock(node, "matlab")) {
                    matlabBlock = node;
                    break;
                }
                if (isBoundary(node) || !isGeneratedJuliaOutput(node)) {
                    return;
                }
                generatedOutput.push(node);
                node = node.nextElementSibling;
            }

            if (!matlabBlock) {
                return;
            }

            var wrapper = document.createElement("div");
            wrapper.className = "code-comparison";
            var juliaColumn = document.createElement("div");
            juliaColumn.className = "code-column code-column-julia";
            var matlabColumn = document.createElement("div");
            matlabColumn.className = "code-column code-column-matlab";

            juliaBlock.parentNode.insertBefore(wrapper, juliaBlock);
            wrapper.appendChild(juliaColumn);
            wrapper.appendChild(matlabColumn);
            juliaColumn.appendChild(juliaBlock);
            generatedOutput.forEach(function (output) {
                juliaColumn.appendChild(output);
            });
            matlabColumn.appendChild(matlabBlock);
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", pairExamples);
    } else {
        pairExamples();
    }
})();
