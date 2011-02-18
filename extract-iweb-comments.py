#!/usr/bin/python
# $Id$
# extract-iweb-comments.py - Dumps iWeb blog comments into individual HTML files for use as HTML snippet code
# Copyright (C) 2009 Benjamin Burke - http://bburke.galvanist.us/
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import os
import gzip
import xml.dom.minidom
import htmlentitydefs
import time
import datetime

domainPath = os.path.expanduser('~/Library/Application Support/iWeb/Domain.sites2')

########################################
#### Functions #########################
########################################

def unicodeToEntities(inputString):
    resultString = u""
    for inputChar in inputString:
         if ord(inputChar) in htmlentitydefs.codepoint2name:
             name = htmlentitydefs.codepoint2name.get(ord(inputChar))
             resultString += "&" + name + ";"
         else:
             resultString += inputChar
    return resultString


def printComment(commentAuthor, commentDate, commentBody):
    return """\
<div id="widget1-item" style="position: relative">
    <div style="padding: 20px 0 0px 0;">
        <div id="widget1-author" class="Comment_Author">%(commentAuthor)s</div>
    </div>

    <div style="padding: 10px 0 10px 0;">
        <div id="widget1-body" class="Comment_Body">%(commentBody)s</div>
    </div>

    <div id="widget1-date" class="Comment_Posted_Date">%(commentDate)s</div>
</div>

<div id="widget1-separator" class="Comment_Separator">
    <div style="border-top: 1px solid #ccc; padding: 0px 0 10px 0;"></div>
</div>
""" % locals()


########################################
#### Main ##############################
########################################

for path, dirs, files in os.walk(domainPath):
    cwd = os.path.basename(path)
    if cwd[0:10] == 'site-blog-':
        for subdir in dirs:
            if subdir[0:10] == "site-page-":
                blogPageArchive = os.path.join(path, subdir, subdir + ".xml.gz")
                blogPage = gzip.open(blogPageArchive);
                dom = xml.dom.minidom.parse(blogPage)
                blogPage.close()
                
                page = dom.getElementsByTagName('bl:site-page')
                pageName = page[0].attributes['sf:name'].value
                headerPrinted = 0
                
                outputFile = open(os.path.expanduser("~/Desktop/%s.html" % (pageName)), 'w')

                outputFile.write("<div id=\"widget1-content\">")
                outputFile.write("<div id=\"widget1-header\" class=\"Comment_Header\" style=\"display: inline; \">")
                outputFile.write("<div style=\"border-bottom: 5px solid #ccc; padding: 0px 0 10px 0;\">")
                outputFile.write("<span class=\"comment-value-comment-count\">COMMENTS ARCHIVED</span></div></div>")

                for comment in dom.getElementsByTagName('sfa:comment'):
                    commentAuthor = comment.attributes['bl:comment-author'].value
                    commentDate = comment.attributes['bl:comment-comparison-date'].value
                    commentBody = comment.childNodes[0].attributes['sfa:string'].value
                    commentDateTime = datetime.datetime.utcfromtimestamp(float(commentDate) + 978307200.0)
                    outputFile.write(printComment(unicodeToEntities(commentAuthor), unicodeToEntities(commentDateTime.strftime("%A, %B %d, %Y - %I:%M %p")), commentBody))
                outputFile.write("</div>")
                outputFile.close()          
          
          