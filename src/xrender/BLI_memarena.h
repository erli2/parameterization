/* 
 * Memory arena ADT
 * 
 * $Id: BLI_memarena.h,v 1.6 2005/03/28 21:49:48 zuster Exp $
 *
 * ***** BEGIN GPL/BL DUAL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. The Blender
 * Foundation also sells licenses for use in proprietary software under
 * the Blender License.  See http://www.blender.org/BL/ for information
 * about this.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL/BL DUAL LICENSE BLOCK *****
 * 
 * Memory arena's are commonly used when the program
 * needs to quickly allocate lots of little bits of
 * data, which are all freed at the same moment.
 * 
 */

#ifndef BLI_MEMARENA_H
#define BLI_MEMARENA_H
#include "DNA_listbase.h"
	/* A reasonable standard buffer size, big
	 * enough to not cause much internal fragmentation, 
	 * small enough not to waste resources
	 */
#define BLI_MEMARENA_STD_BUFSIZE	(1<<14)

struct MemArena;
typedef struct MemArena MemArena;


struct MemArena*	BLI_memarena_new	(int bufsize);
void				BLI_memarena_free	(struct MemArena *ma);

void*				BLI_memarena_alloc	(struct MemArena *ma, int size);

//struct MemArena;

typedef void (*LinkNodeFreeFP)(void *link);
typedef void (*LinkNodeApplyFP)(void *link);

struct LinkNode;
typedef struct LinkNode {
	struct LinkNode *next;
	void *link;
} LinkNode;

int		BLI_linklist_length		(struct LinkNode *list);

void	BLI_linklist_reverse	(struct LinkNode **listp);

void	BLI_linklist_prepend		(struct LinkNode **listp, void *ptr);
void	BLI_linklist_append	    	(struct LinkNode **listp, void *ptr);
void	BLI_linklist_prepend_arena	(struct LinkNode **listp, void *ptr, struct MemArena *ma);

void	BLI_linklist_free		(struct LinkNode *list, LinkNodeFreeFP freefunc);
void	BLI_linklist_apply		(struct LinkNode *list, LinkNodeApplyFP applyfunc);

void	BLI_addtail(ListBase *listbase, void *vlink);
void	BLI_freelistN(ListBase *listbase);

#endif

