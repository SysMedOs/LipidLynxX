from flask import abort
from flask_restful import Resource, fields, marshal_with

from .parsers import convert_get_parser, convert_post_parser


nested_tag_fields = {
    'id': fields.Integer(),
    'title': fields.String()
}

post_fields = {
    'id': fields.Integer(),
    'title': fields.String(),
}


class ConvertApi(Resource):
    @marshal_with(post_fields)
    def get(self, lynx_id=None):
        if lynx_id:
            post = Post.query.get(lynx_id)
            if not post:
                abort(404)

            return post
        else:
            args = convert_get_parser.parse_args()
            lynx_id = args['id'] or ''

            if lynx_id:
                pass
            else:
                pass

            return posts.items

    def post(self, post_id=None):
        if post_id:
            abort(400)
        else:
            args = post_post_parser.parse_args(strict=True)

            user = User.verify_auth_token(args['token'])
            if not user:
                abort(401)

            new_post = Post(args['title'])
            new_post.user = user
            new_post.date = datetime.datetime.now()
            new_post.text = args['text']

            if args['tags']:
                for item in args['tags']:
                    tag = Tag.query.filter_by(title=item).first()

                    # Add the tag if it exists. If not, make a new tag
                    if tag:
                        new_post.tags.append(tag)
                    else:
                        new_tag = Tag(item)
                        new_post.tags.append(new_tag)

            db.session.add(new_post)
            db.session.commit()
            return new_post.id, 201
