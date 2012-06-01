var NAVTREE =
[
  [ "frodo-l2-pipeline", "index.html", [
    [ "File List", "files.html", [
      [ "FITS2jpeg.c", "_f_i_t_s2jpeg_8c.html", null ],
      [ "frodo_aux_numexts.c", "frodo__aux__numexts_8c.html", null ],
      [ "frodo_aux_numexts.h", "frodo__aux__numexts_8h.html", null ],
      [ "frodo_aux_spec2tsv.c", "frodo__aux__spec2tsv_8c.html", null ],
      [ "frodo_aux_spec2tsv.h", "frodo__aux__spec2tsv_8h.html", null ],
      [ "frodo_config.h", "frodo__config_8h.html", null ],
      [ "frodo_error_handling.c", "frodo__error__handling_8c.html", null ],
      [ "frodo_error_handling.h", "frodo__error__handling_8h.html", null ],
      [ "frodo_functions.c", "frodo__functions_8c.html", null ],
      [ "frodo_functions.h", "frodo__functions_8h.html", null ],
      [ "frodo_red_arcfit.c", "frodo__red__arcfit_8c.html", null ],
      [ "frodo_red_arcfit.h", "frodo__red__arcfit_8h.html", null ],
      [ "frodo_red_arcfit2.c", "frodo__red__arcfit2_8c.html", null ],
      [ "frodo_red_arcfit2.h", "frodo__red__arcfit2_8h.html", null ],
      [ "frodo_red_correct_throughput.c", "frodo__red__correct__throughput_8c.html", null ],
      [ "frodo_red_correct_throughput.h", "frodo__red__correct__throughput_8h.html", null ],
      [ "frodo_red_extract_simple.c", "frodo__red__extract__simple_8c.html", null ],
      [ "frodo_red_extract_simple.h", "frodo__red__extract__simple_8h.html", null ],
      [ "frodo_red_findpeaks_simple.c", "frodo__red__findpeaks__simple_8c.html", null ],
      [ "frodo_red_findpeaks_simple.h", "frodo__red__findpeaks__simple_8h.html", null ],
      [ "frodo_red_findpeaks_simple_clean.c", "frodo__red__findpeaks__simple__clean_8c.html", null ],
      [ "frodo_red_findpeaks_simple_clean.h", "frodo__red__findpeaks__simple__clean_8h.html", null ],
      [ "frodo_red_rebin.c", "frodo__red__rebin_8c.html", null ],
      [ "frodo_red_rebin.h", "frodo__red__rebin_8h.html", null ],
      [ "frodo_red_reformat.c", "frodo__red__reformat_8c.html", null ],
      [ "frodo_red_reformat.h", "frodo__red__reformat_8h.html", null ],
      [ "frodo_red_subsky.c", "frodo__red__subsky_8c.html", null ],
      [ "frodo_red_subsky.h", "frodo__red__subsky_8h.html", null ],
      [ "frodo_red_trace.c", "frodo__red__trace_8c.html", null ],
      [ "frodo_red_trace.h", "frodo__red__trace_8h.html", null ],
      [ "imcopy.c", "imcopy_8c.html", null ],
      [ "jpegsubs.c", "jpegsubs_8c.html", null ],
      [ "jpegsubs.h", "jpegsubs_8h.html", null ],
      [ "liststruc.c", "liststruc_8c.html", null ],
      [ "modhead.c", "modhead_8c.html", null ]
    ] ],
    [ "Globals", "globals.html", null ]
  ] ]
];

function createIndent(o,domNode,node,level)
{
  if (node.parentNode && node.parentNode.parentNode)
  {
    createIndent(o,domNode,node.parentNode,level+1);
  }
  var imgNode = document.createElement("img");
  if (level==0 && node.childrenData)
  {
    node.plus_img = imgNode;
    node.expandToggle = document.createElement("a");
    node.expandToggle.href = "javascript:void(0)";
    node.expandToggle.onclick = function() 
    {
      if (node.expanded) 
      {
        $(node.getChildrenUL()).slideUp("fast");
        if (node.isLast)
        {
          node.plus_img.src = node.relpath+"ftv2plastnode.png";
        }
        else
        {
          node.plus_img.src = node.relpath+"ftv2pnode.png";
        }
        node.expanded = false;
      } 
      else 
      {
        expandNode(o, node, false);
      }
    }
    node.expandToggle.appendChild(imgNode);
    domNode.appendChild(node.expandToggle);
  }
  else
  {
    domNode.appendChild(imgNode);
  }
  if (level==0)
  {
    if (node.isLast)
    {
      if (node.childrenData)
      {
        imgNode.src = node.relpath+"ftv2plastnode.png";
      }
      else
      {
        imgNode.src = node.relpath+"ftv2lastnode.png";
        domNode.appendChild(imgNode);
      }
    }
    else
    {
      if (node.childrenData)
      {
        imgNode.src = node.relpath+"ftv2pnode.png";
      }
      else
      {
        imgNode.src = node.relpath+"ftv2node.png";
        domNode.appendChild(imgNode);
      }
    }
  }
  else
  {
    if (node.isLast)
    {
      imgNode.src = node.relpath+"ftv2blank.png";
    }
    else
    {
      imgNode.src = node.relpath+"ftv2vertline.png";
    }
  }
  imgNode.border = "0";
}

function newNode(o, po, text, link, childrenData, lastNode)
{
  var node = new Object();
  node.children = Array();
  node.childrenData = childrenData;
  node.depth = po.depth + 1;
  node.relpath = po.relpath;
  node.isLast = lastNode;

  node.li = document.createElement("li");
  po.getChildrenUL().appendChild(node.li);
  node.parentNode = po;

  node.itemDiv = document.createElement("div");
  node.itemDiv.className = "item";

  node.labelSpan = document.createElement("span");
  node.labelSpan.className = "label";

  createIndent(o,node.itemDiv,node,0);
  node.itemDiv.appendChild(node.labelSpan);
  node.li.appendChild(node.itemDiv);

  var a = document.createElement("a");
  node.labelSpan.appendChild(a);
  node.label = document.createTextNode(text);
  a.appendChild(node.label);
  if (link) 
  {
    a.href = node.relpath+link;
  } 
  else 
  {
    if (childrenData != null) 
    {
      a.className = "nolink";
      a.href = "javascript:void(0)";
      a.onclick = node.expandToggle.onclick;
      node.expanded = false;
    }
  }

  node.childrenUL = null;
  node.getChildrenUL = function() 
  {
    if (!node.childrenUL) 
    {
      node.childrenUL = document.createElement("ul");
      node.childrenUL.className = "children_ul";
      node.childrenUL.style.display = "none";
      node.li.appendChild(node.childrenUL);
    }
    return node.childrenUL;
  };

  return node;
}

function showRoot()
{
  var headerHeight = $("#top").height();
  var footerHeight = $("#nav-path").height();
  var windowHeight = $(window).height() - headerHeight - footerHeight;
  navtree.scrollTo('#selected',0,{offset:-windowHeight/2});
}

function expandNode(o, node, imm)
{
  if (node.childrenData && !node.expanded) 
  {
    if (!node.childrenVisited) 
    {
      getNode(o, node);
    }
    if (imm)
    {
      $(node.getChildrenUL()).show();
    } 
    else 
    {
      $(node.getChildrenUL()).slideDown("fast",showRoot);
    }
    if (node.isLast)
    {
      node.plus_img.src = node.relpath+"ftv2mlastnode.png";
    }
    else
    {
      node.plus_img.src = node.relpath+"ftv2mnode.png";
    }
    node.expanded = true;
  }
}

function getNode(o, po)
{
  po.childrenVisited = true;
  var l = po.childrenData.length-1;
  for (var i in po.childrenData) 
  {
    var nodeData = po.childrenData[i];
    po.children[i] = newNode(o, po, nodeData[0], nodeData[1], nodeData[2],
        i==l);
  }
}

function findNavTreePage(url, data)
{
  var nodes = data;
  var result = null;
  for (var i in nodes) 
  {
    var d = nodes[i];
    if (d[1] == url) 
    {
      return new Array(i);
    }
    else if (d[2] != null) // array of children
    {
      result = findNavTreePage(url, d[2]);
      if (result != null) 
      {
        return (new Array(i).concat(result));
      }
    }
  }
  return null;
}

function initNavTree(toroot,relpath)
{
  var o = new Object();
  o.toroot = toroot;
  o.node = new Object();
  o.node.li = document.getElementById("nav-tree-contents");
  o.node.childrenData = NAVTREE;
  o.node.children = new Array();
  o.node.childrenUL = document.createElement("ul");
  o.node.getChildrenUL = function() { return o.node.childrenUL; };
  o.node.li.appendChild(o.node.childrenUL);
  o.node.depth = 0;
  o.node.relpath = relpath;

  getNode(o, o.node);

  o.breadcrumbs = findNavTreePage(toroot, NAVTREE);
  if (o.breadcrumbs == null)
  {
    o.breadcrumbs = findNavTreePage("index.html",NAVTREE);
  }
  if (o.breadcrumbs != null && o.breadcrumbs.length>0)
  {
    var p = o.node;
    for (var i in o.breadcrumbs) 
    {
      var j = o.breadcrumbs[i];
      p = p.children[j];
      expandNode(o,p,true);
    }
    p.itemDiv.className = p.itemDiv.className + " selected";
    p.itemDiv.id = "selected";
    $(window).load(showRoot);
  }
}

