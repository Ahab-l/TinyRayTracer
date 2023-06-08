
#include <fstream>


#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hitablelist.h"
#include "sphere.h"
#include "material.h"
#include <iostream>
#include"moving_sphere.h"
#include"bvh.h"
#include"aarect.h"
#include"box.h"
#include "constant_medium.h"
#include"pdf.h"
hittable_list cornell_smoke() {
	hittable_list objects;

	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(7, 7, 7));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	//objects.add(make_shared<zx_rect>(227, 332, 213, 343, 554, light));
	objects.add(make_shared<zx_rect>(127, 432, 113, 443, 554, light));
	objects.add(make_shared<zx_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<zx_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, vec3(265, 0, 295));

	shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, vec3(130, 0, 65));

	objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
	objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));

	return objects;
}
hittable_list cornell_box_only_lamb() {
	hittable_list objects;

	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(15, 15, 15));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));

	objects.add(make_shared<flip_face>(make_shared<zx_rect>(227, 332, 213, 343, 554, light)));
	objects.add(make_shared<zx_rect>(0, 555, 0, 555, 0, white));
	//objects.add(make_shared<zx_rect>(0, 555, 0, 555, 278, white));
	objects.add(make_shared<zx_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

	shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), aluminum);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, vec3(265, 0, 295));
	objects.add(box1);

	//shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
	//box2 = make_shared<rotate_y>(box2, -18);
	//box2 = make_shared<translate>(box2, vec3(130, 0, 65));
	//objects.add(box2);
	auto glass = make_shared<dielectric>(1.5);
	objects.add(make_shared<sphere>(point3(190, 90, 190), 90, glass));

	return objects;
}
hittable_list cornell_box() {
	hittable_list objects;

	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(15, 15, 15));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<zx_rect>(227, 332, 213, 343,  554, light));
	objects.add(make_shared<zx_rect>( 0, 555, 0, 555, 0,  white));
	objects.add(make_shared<zx_rect>(0, 555, 0, 555,  555, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, vec3(265, 0, 295));
	objects.add(box1);

	shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, vec3(130, 0, 65));
	objects.add(box2);
	return objects;
}
hittable_list two_spheres() {
	hittable_list objects;

	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

	objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));
	objects.add(make_shared<sphere>(point3(0, 10, 0), 10, make_shared<lambertian>(checker)));

	return objects;
}
hittable_list two_perlin_spheres() {
	hittable_list objects;

	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	return objects;
}
hittable_list earth_map_spheres() {
	hittable_list objects;

	auto text_map = make_shared<image_texture>("earthmap.jpg");
	objects.add(make_shared<sphere>(point3(0, 3, 0), 5, make_shared<lambertian>(text_map)));

	return objects;
}
hittable_list simple_light() {
	hittable_list objects;

	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
	objects.add(make_shared<xy_rect>(1, 3, 1, 3, -3, difflight));
	objects.add(make_shared<yz_rect>(1, 3, 1, 3, -4, difflight));
	objects.add(make_shared<zx_rect>(-1, 1, -1, 1, 5, difflight));

	return objects;
}

hittable_list random_scene() {
	hittable_list world;

	
	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
	world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			auto choose_mat = random_double();
			point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

			if ((center - point3(4, 0.2, 0)).length() > 0.9) {
				shared_ptr<material> sphere_material;


				if (choose_mat < 0.8) {
					// diffuse
					auto albedo = random() * random();
					sphere_material = make_shared<lambertian>(albedo);
					auto center2 = center + vec3(0, random_double(0, .5), 0);
					world.add(make_shared<moving_sphere>(
						center, center2, 0.0, 1.0, 0.2, sphere_material));
				}
				else if (choose_mat < 0.95) {
					// metal
					auto albedo = random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<sphere>(center, 0.2, sphere_material));
				}
				else {
					// glass
					sphere_material = make_shared<dielectric>(1.5);
					world.add(make_shared<sphere>(center, 0.2, sphere_material));
				}
			}
		}
	}
	auto material_emit = make_shared<diffuse_light>(color(1.0, 1.0, 1.0));



	auto material1 = make_shared<dielectric>(1.5);
	world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material_emit));

	auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
	world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

	auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
	world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

	return world;
}
color ray_color_bvh(const ray& r,const hittable& world,int depth,const bvh_node &bn) {
	hit_record rec;
	scatter_record srec;
	double pdf;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}


	if (bn.hit(r, 0.00001, infinity, rec)) {
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, srec))
			return srec.attenuation * ray_color_bvh(srec.specular_ray, world, depth - 1,bn);
		return color(0, 0, 0);

	}

	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}
color ray_color_bvh_emit(const ray& r, const color& background,const hittable& world, int depth, const bvh_node& bn) {
	hit_record rec;
	scatter_record srec;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}

	if (!bn.hit(r, 0.001, infinity, rec))
		return background;


	color attenuation;
	color emitted = rec.mat_ptr->emitted(r,rec,rec.u, rec.v, rec.p);

	if (!rec.mat_ptr->scatter(r,rec,srec))
	{
		return emitted;
	}

	//return emitted + albedo * ray_color_bvh_emit(scattered, background, world, depth - 1,bn);
	return emitted
		+ srec.attenuation *rec.mat_ptr->scattering_pdf(r, rec, srec.specular_ray) * ray_color_bvh_emit(srec.specular_ray, background, world, depth - 1,bn) / (dot(rec.normal, srec.specular_ray.direction()) / pi);

}color ray_color_light_sample(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn) {
	if (depth == 48)
		std::cout << 1;
	hit_record rec;
	scatter_record srec;
	double pdf;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}

	if (!bn.hit(r, 0.001, infinity, rec))
		return background;

	color attenuation;
	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);

	if (!rec.mat_ptr->scatter(r, rec, srec))
	{
		return emitted;
	}

	auto on_light = point3(random_double(213, 343), 554, random_double(227, 332));
	auto to_light = on_light - rec.p;
	auto distance_squared = to_light.length_squared();
	to_light = unit_vector(to_light);

	if (dot(to_light, rec.normal) < 0)
		return emitted;

	double light_area = (343 - 213) * (332 - 227);
	auto light_cosine = fabs(to_light.y());
	if (light_cosine < 0.000001)
		return emitted;

	pdf = distance_squared / (light_cosine * light_area);
	srec.specular_ray = ray(rec.p, to_light, r.time());

	//return emitted + albedo * ray_color_bvh_emit(scattered, background, world, depth - 1,bn);
	return emitted
		+ srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, srec.specular_ray)* ray_color_light_sample(srec.specular_ray, background, world, depth - 1, bn) / pdf;

}

color ray_color_cosine_pdf(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn) {
	hit_record rec;
	double pdf;

	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}

	if (!bn.hit(r, 0.001, infinity, rec))
		return background;
	scatter_record srec;
	ray scattered;
	color attenuation;
	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
	color albedo;
	double pdf_val;
	if (!rec.mat_ptr->scatter(r, rec, srec))
	{
		return emitted;
	}
	cosine_pdf p(rec.normal);
	scattered = ray(rec.p, p.generate(), r.time());
	pdf_val = p.value(scattered.direction());
	
	double x = rec.mat_ptr->scattering_pdf(r, rec, scattered);
	//return emitted + albedo * ray_color_bvh_emit(scattered, background, world, depth - 1,bn);
	return emitted
		+ albedo*x
				 * ray_color_cosine_pdf(scattered, background, world, depth - 1, bn) /pdf_val;

}
//color ray_color_mixed(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn,shared_ptr<hittable>&lights) {
//	hit_record rec;
//	double pdf;
//
//	// If we've exceeded the ray bounce limit, no more light is gathered.
//	if (depth <= 0) {
//		return color(0, 0, 0);
//	}
//
//	if (!bn.hit(r, 0.001, infinity, rec))
//		return background;
//
//	ray scattered;
//	color attenuation;
//	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
//	color albedo;
//	double pdf_val;
//
//	if (!rec.mat_ptr->scatter(r, rec, albedo, scattered, pdf_val))
//	{
//		return emitted;
//	}
//
//	auto p0 = make_shared<hittable_pdf>(lights,rec.p);
//	auto p1 = make_shared<cosine_pdf>(rec.normal);
//	mixture_pdf mixed_pdf(p0, p1);
//
//	scattered = ray(rec.p, mixed_pdf.generate(), r.time());
//	pdf_val = mixed_pdf.value(scattered.direction());
//
//	double x = rec.mat_ptr->scattering_pdf(r, rec, scattered);
//	//return emitted + albedo * ray_color_bvh_emit(scattered, background, world, depth - 1,bn);
//	return emitted
//		+ albedo * x
//		* ray_color_mixed(scattered, background, world, depth - 1, bn,lights) / pdf_val;
//
//}
//color ray_color_mixed(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn, shared_ptr<hittable>& lights) {
//	hit_record rec;
//	double pdf;
//
//	// If we've exceeded the ray bounce limit, no more light is gathered.
//	if (depth <= 0) {
//		return color(0, 0, 0);
//	}
//
//	if (!bn.hit(r, 0.001, infinity, rec))
//		return background;
//
//	scatter_record srec;
//	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
//	if (!rec.mat_ptr->scatter(r, rec, srec))
//		return emitted;
//	if (srec.is_specular) {
//		return srec.attenuation * ray_color_mixed(srec.specular_ray, background, world,depth-1, bn, lights);
//	}
//	auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
//
//	mixture_pdf p(light_ptr, srec.pdf_ptr);
//
//	ray scattered = ray(rec.p, p.generate(), r.time());
//	auto pdf_val = p.value(scattered.direction());
//
//	return emitted + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
//		* ray_color_mixed(scattered, background, world, depth - 1,bn,lights) / pdf_val;
//
//
//}
color ray_color_mixed(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn, shared_ptr<hittable> lights) {
	hit_record rec;
	double pdf;

	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}

	if (!bn.hit(r, 0.001, infinity, rec))
		return background;

	scatter_record srec;
	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
	if (!rec.mat_ptr->scatter(r, rec, srec))
		return emitted;
	if (srec.is_specular) {
		return srec.attenuation * ray_color_mixed(srec.specular_ray, background, world, depth - 1, bn, lights);
	}
	//auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);

	//mixture_pdf p(light_ptr, srec.pdf_ptr);
	auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
	mixture_pdf p(light_ptr, srec.pdf_ptr);

	ray scattered = ray(rec.p, p.generate(), r.time());
	auto pdf_val = p.value(scattered.direction());

	return emitted + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
		* ray_color_mixed(scattered, background, world, depth - 1, bn, lights) / pdf_val;


}
color ray_color_hitt(const ray& r, const color& background, const hittable& world, int depth, const bvh_node& bn,shared_ptr<hittable>&lights) {
	hit_record rec;
	scatter_record srec;
	double pdf;

	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}

	if (!bn.hit(r, 0.001, infinity, rec))
		return background;

	ray scattered;
	color attenuation;
	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
	color albedo;
	double pdf_val;
	if (!rec.mat_ptr->scatter(r, rec, srec))
	{
		return emitted;
	}
	hittable_pdf light_pdf(lights, rec.p);
	srec.specular_ray = ray(rec.p, light_pdf.generate(), r.time());
	pdf_val = light_pdf.value(srec.specular_ray.direction());
	
	//double x = rec.mat_ptr->scattering_pdf(r, rec, scattered);
	//return emitted + albedo * ray_color_bvh_emit(scattered, background, world, depth - 1,bn);
	return emitted
		+ srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, srec.specular_ray)
		* ray_color_hitt(srec.specular_ray, background, world, depth - 1, bn,lights) / pdf_val;

}

color ray_color(const ray& r, const hittable& world, int depth, const bvh_node& bn) {
	hit_record rec;
	scatter_record srec;
	double pdf;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) {
		return color(0, 0, 0);
	}


	if (world.hit(r, 0.00001, infinity, rec)) {
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec,srec))
			return srec.attenuation * ray_color(srec.specular_ray, world, depth - 1, bn);
		return color(0, 0, 0);

	}

	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

int main() {


	// Image

	auto aspect_ratio = 21.0 / 9.0;
	int image_width = 100;
	int samples_per_pixel = 1000;
	const int max_depth         = 50;
	// World
	hittable_list world;
	//objects.add(make_shared<flip_face>(make_shared<zx_rect>(227, 332, 213, 343, 554, light)));
	auto lights = make_shared<hittable_list>();
	lights->add(make_shared<zx_rect>(227, 332, 213, 343, 554, shared_ptr<material>()));
	//lights.add(make_shared<sphere>(point3(190, 90, 190), 90, shared_ptr<material>()));
	//shared_ptr<hittable> lights = //make_shared<zx_rect>(227, 332, 213, 343, 554, shared_ptr<material>());
								  //make_shared<sphere>(point3(190, 90, 190), 90, shared_ptr<material>());
	//auto lights = make_shared<hittable_list>();
	//lights->add(make_shared<zx_rect>(227, 332, 213, 343, 554, shared_ptr<material>()));
	//lights->add(make_shared<sphere>(point3(190, 90, 190), 90, shared_ptr<material>()));
	point3 lookfrom;
	point3 lookat;
	auto vfov = 40.0;
	auto aperture = 0.0;
	color background(0, 0, 0);

	auto time0 = 0.0;
	auto time1 = 1.0;

	switch (9) {
	case 1:
		world = random_scene();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		aperture = 0.1;
		break;
	case 2:
		world = two_spheres();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;
	case 3:
		world = two_perlin_spheres();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;

	case 4:
		world = earth_map_spheres();
		lookfrom = point3(0, 2, 10);
		lookat = point3(0, 2, 0);
		vfov = 90.0;
		break;

	case 5:
		background = color(0.0, 0.0, 0.0);
		world = random_scene();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;
	case 6:
		world = simple_light();
		samples_per_pixel = 400;
		background = color(0, 0, 0);
		lookfrom = point3(26, 3, 6);
		lookat = point3(0, 2, 0);
		vfov = 20.0;
		break;
	case 7:
		world = cornell_box();
		aspect_ratio = 1.0;
		image_width = 200;
		samples_per_pixel = 200;
		background = color(0, 0, 0);
		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		break;
	case 8:
		world = cornell_smoke();
		aspect_ratio = 1.0;
		image_width = 600;
		samples_per_pixel = 200;
		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		break;
	case 9: 
		world = cornell_box_only_lamb();
		image_width = 300;
		aspect_ratio = 1.0;
		samples_per_pixel = 200;
		

		color background(0, 0, 0);

		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;

		time0 = 0.0;
		time1 = 1.0;
		break;

	}


	bvh_node bn(world, 0.0, 1.0);
	// Camera

	vec3 vup(0, 1, 0);
	auto dist_to_focus = 10.0;
	int image_height = static_cast<int>(image_width / aspect_ratio);

	camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, time0, time1);

	
	// Render

	std::ofstream outimage;
	//outimage.open("cornel_pdf_mixed_result.ppm");
	outimage.open("cornel_pdf_samplelight1.ppm");
	outimage << "P3\n" << image_width << " " << image_height << "\n255\n";

	#pragma omp parallel
	#pragma omp for
	for (int j = image_height - 1; j >= 0; --j) {
		color pixel_color(0, 0, 0);
		for (int i = 0; i < image_width; ++i) {
			color pixel_color(0, 0, 0);
			for ( int s = 0; s < samples_per_pixel; ++s)
			{
				auto u = double(i+random_double()) / (image_width - 1);
				auto v = double(j+random_double()) / (image_height - 1);
				ray  r = cam.get_ray(u, v);
				pixel_color += ray_color_mixed(r, background,world,max_depth,bn,lights);
			}
			write_color(outimage, pixel_color,samples_per_pixel);
			std::cout << "pixel:in line" << j << "row" << i << "has been computed"<<std::endl;
		}
	}

	std::cerr << "\nDone.\n";
}