/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xml.h
 * @brief  declarations for XML parsing
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_XML_H__
#define __SCIP_XML_H__

#ifdef __cplusplus
extern "C" {
#endif


typedef struct XML_ATTR_struct XML_ATTR;

struct XML_ATTR_struct
{
   char*                 name;
   char*                 value;
   XML_ATTR*             next;
};

typedef struct XML_NODE_struct XML_NODE;

struct XML_NODE_struct
{
   char*                 name;
   int                   lineno;
   XML_ATTR*             attrlist;
   XML_NODE*             parent;
   XML_NODE*             prevsibl;
   XML_NODE*             nextsibl;
   XML_NODE*             firstchild;
   XML_NODE*             lastchild;
   char*                 data;           /* does not come together with children */
};

/** Parse file */
extern
XML_NODE* xmlProcess(
   const char*           filename            /**< XML file name */
   );

/** create new node */
extern
XML_NODE* xmlNewNode(
   const char*           name,
   int                   lineno
   );

/** create new attribute */
extern
XML_ATTR* xmlNewAttr(
   const char*           name,
   const char*           value
   );

/** add attribute */
extern
void xmlAddAttr(
   XML_NODE*             n,
   XML_ATTR*             a
   );

/** append child node */
extern
void xmlAppendChild(
   XML_NODE*             parent,
   XML_NODE*             child
   );

/** free node */
extern
void xmlFreeNode(
   XML_NODE*             node
   );

/** output node */
extern
void xmlShowNode(
   const XML_NODE*       root
   );

/** get attribute value */
extern
const char* xmlGetAttrval(
   const XML_NODE*       node,
   const char*           name
   );

/** return first node */
extern
const XML_NODE* xmlFirstNode(
   const XML_NODE*       node,
   const char*           name
   );

/** return next node */
extern
const XML_NODE* xmlNextNode(
   const XML_NODE*       node,
   const char*           name
   );

/** find node */
extern
const XML_NODE* xmlFindNode(
   const XML_NODE*       node,
   const char*           name
   );

/** find node with bound on the depth */
extern
const XML_NODE* xmlFindNodeMaxdepth(
   const XML_NODE*       node,               /**< current node - use start node to begin */
   const char*           name,               /**< name of tag to search for */
   int                   depth,              /**< current depth - start with 0 */
   int                   maxdepth            /**< maximal depth */
   );

/** return next sibling */
extern
const XML_NODE* xmlNextSibl(
   const XML_NODE*       node
   );

/** return previous sibling */
extern
const XML_NODE* xmlPrevSibl(
   const XML_NODE*       node
   );

/** return first child */
extern
const XML_NODE* xmlFirstChild(
   const XML_NODE*       node
   );

/** return last child */
extern
const XML_NODE* xmlLastChild(
   const XML_NODE*       node
   );

/** return name of node */
extern
const char* xmlGetName(
   const XML_NODE*       node
   );

/** get line number */
extern
int xmlGetLine(
   const XML_NODE*       node
   );

/** get data */
extern
const char* xmlGetData(
   const XML_NODE*       node
   );

/** find PCDATA */
extern
const char* xmlFindPcdata(
   const XML_NODE*       node,
   const char*           name
   );

#ifdef __cplusplus
}
#endif

#endif
