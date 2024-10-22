/// Modified version of ROOT-version selector
/// Modified by Bernd Riederer
/// (c) Axel Naumann, CERN; 2020-03-02

/// Configurable section.

// What the master is called. Leave untouched if master has no doc.
let master = 'master';

/// Pathname part of the URL containing the different versioned doxygen
/// subdirectories. Must be browsable.
let urlroot = '/';

///=============================================================================
// Remove trailing '/'.
if (urlroot[urlroot.length - 1] == '/')
  urlroot = urlroot.substring(0, urlroot.length - 1)
let urlrootdirs = urlroot.split('/').length;

function url2version(patharr) {
  // Given the directory array of a URL (i.e. without domain), extract the
  // version corresponding to the URL.
  // E.g. for `https://example.com/doc/master/classX.html`, the directory array
  // becomes `["doc", "master, "classX.html"]. This function might return
  // the second element, `master`.
  return patharr[urlrootdirs];
}

let patharr = window.location.pathname.replace(/\/+/g, '/').split('/');
let thisvers = url2version(patharr);
$('.dropbtn').html(thisvers);

// https://stackoverflow.com/questions/30622369
$.get(urlroot + '/', (data) => {
  let ret = parseDirectoryListing(data);
  $('.dropdown-content').append(ret.join(''));

  // Now check which links actually exist, and remove the href for those
  // that don't.
  $('.dropdown-content')
    .find('.verslink')
    .each(function () {
      var el = $(this);
      var request = new XMLHttpRequest();
      request.open('HEAD', el.attr('href'), true);
      request.onreadystatechange = function () {
        if (request.readyState === 4) {
          if (request.status === 404) {
            el.removeAttr('href');
            el.css({
              'color': "gray",
              'text-decoration': 'line-through'
            });
          }
        }
      };
      request.send();
    });
});

function parseDirectoryListing(text) {
  let docs = text.match(/href="([^/][^"]+)"/g); // only directories.
  docs = docs.map((x) =>
    x.replace(/\.\//g, '')
      .replace(/href="/, '')
      .replace(/"$/, ''))
    .sort().reverse();
  docs = docs.filter(function (line) {
    // suppress link to current version
    return !RegExp(thisvers).test(line)
      // suppress hidden dirs
      && !/^[.]/.test(line);
  });

  if (docs.includes(master)) {
    // Move "master" first.
    docs.splice(docs.findIndex(el => el == master), 1);
    docs.splice(0, 0, master);
  }
  docs = docs.map((x) => '<a class="verslink" href="'
    + patharr.slice(0, urlrootdirs).join('/')
    + '/' + x + '/' + patharr[patharr.length - 1] + '">'
    + x
    + '</a>');
  return docs;
}
